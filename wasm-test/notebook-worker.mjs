/**
 * PythonSCAD notebook kernel worker: owns Emscripten module and Python REPL.
 */
import factory from './pythonscad.js';
import {MSG, REPL_PRELUDE, buildExportScript, combinePythonOutput, exportPathForFormat, wrapNotebookExec,} from './notebook-protocol.mjs';

let mod = null;
let evalPython = null;
let pythonInited = false;
let jobChain = Promise.resolve();

function errMsg(e)
{
  if (e == null) return 'unknown error';
  if (typeof e === 'string') return e;
  if (e instanceof Error) return e.message || String(e);
  try {
    return String(e.message ?? e);
  } catch (_) {
    return String(e);
  }
}

function enqueueJob(fn)
{
  const run = jobChain.then(fn, fn);
  jobChain = run.catch(() => {});
  return run;
}

function reply(data, transfer)
{
  self.postMessage(data, transfer || []);
}

function replyError(id, e, crashed = false)
{
  reply({type: MSG.ERROR, id, message: errMsg(e), crashed});
}

function makeEvalContext()
{
  const stdoutBuf = [];
  const stderrBuf = [];
  return {
    stdoutBuf,
    stderrBuf,
    clear() {
      stdoutBuf.length = 0;
      stderrBuf.length = 0;
    },
    evaluate(source) {
      this.clear();
      const returnValue = evalPython.evaluate(source, false);
      return combinePythonOutput(stdoutBuf, stderrBuf, returnValue);
    },
  };
}

let evalCtx = null;

async function initKernel(wasmBinary, dataBinary)
{
  mod = null;
  evalPython = null;
  pythonInited = false;
  evalCtx = null;

  const factoryOpts = {
    noInitialRun: true,
    print: (t) => {
      if (evalCtx) evalCtx.stdoutBuf.push(t);
    },
    printErr: (t) => {
      if (evalCtx) evalCtx.stderrBuf.push(t);
    },
  };
  if (wasmBinary) factoryOpts.wasmBinary = wasmBinary;
  if (dataBinary) factoryOpts.getPreloadedPackage = () => dataBinary;

  mod = await factory(factoryOpts);
  evalPython = {
    init: mod.cwrap('EmsInitPython', null, []),
    evaluate: mod.cwrap('EmsEvaluatePython', 'string', ['string', 'number']),
    finish: mod.cwrap('EmsFinishPython', null, []),
  };

  evalPython.init();
  pythonInited = true;
  evalCtx = makeEvalContext();

  evalCtx.evaluate(REPL_PRELUDE);
  evalCtx.clear();
}

async function handleInit(msg)
{
  try {
    await initKernel(msg.wasmBinary, msg.dataBinary);
    reply({type: MSG.RESULT, id: msg.id, ready: true});
  } catch (e) {
    replyError(msg.id, e, false);
  }
}

async function handleEvaluate(msg)
{
  if (!pythonInited || !evalCtx) {
    replyError(msg.id, 'Kernel not initialized');
    return;
  }
  try {
    const output = evalCtx.evaluate(wrapNotebookExec(msg.code));
    reply({type: MSG.RESULT, id: msg.id, output});
  } catch (e) {
    replyError(msg.id, e, true);
  }
}

async function handleExport(msg)
{
  if (!pythonInited || !evalCtx || !mod) {
    replyError(msg.id, 'Kernel not initialized');
    return;
  }

  const format = String(msg.format || 'stl');
  const path = exportPathForFormat(format);

  try {
    try {
      mod.FS.mkdir('/tmp');
    } catch (_) { /* exists */
    }
    try {
      mod.FS.unlink(path);
    } catch (_) { /* absent */
    }

    const output = evalCtx.evaluate(buildExportScript(format, path));
    if (output.includes('NO_MODEL_SHOWN_YET')) {
      reply({type: MSG.EXPORT_RESULT, id: msg.id, status: 'no-model'});
      return;
    }

    const fileBytes = mod.FS.readFile(path);
    const buffer =
      fileBytes.buffer.slice(fileBytes.byteOffset, fileBytes.byteOffset + fileBytes.byteLength);
    reply(
      {
        type: MSG.EXPORT_RESULT,
        id: msg.id,
        status: 'ok',
        format,
        byteLength: buffer.byteLength,
        fileBytes: buffer,
      },
      [buffer]);
  } catch (e) {
    replyError(msg.id, e, true);
  }
}

self.onmessage = (event) => {
  const msg = event.data;
  if (!msg || typeof msg.type !== 'string') return;

  if (msg.type === MSG.INIT) {
    enqueueJob(() => handleInit(msg));
    return;
  }
  if (msg.type === MSG.EVALUATE) {
    enqueueJob(() => handleEvaluate(msg));
    return;
  }
  if (msg.type === MSG.EXPORT) {
    enqueueJob(() => handleExport(msg));
  }
};
