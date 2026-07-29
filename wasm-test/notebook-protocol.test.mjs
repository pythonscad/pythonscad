import assert from 'node:assert/strict';
import test from 'node:test';

import {MSG, REPL_PRELUDE, autoShowLastLine, buildExportScript, combinePythonOutput, exportPathForFormat, resolveWorkerResponse, wrapNotebookExec,} from './notebook-protocol.mjs';
import {NotebookKernel} from './notebook-kernel.mjs';

test('REPL prelude includes notebook hooks', () => {
  assert.match(REPL_PRELUDE, /PRELUDE_LOADED_OK/);
  assert.match(REPL_PRELUDE, /_notebook_exec/);
});

test('autoShowLastLine wraps bare identifiers', () => {
  assert.equal(autoShowLastLine('mymodel\n'), 'show(mymodel)\n');
  assert.equal(autoShowLastLine('pass\n'), 'pass\n');
  assert.equal(autoShowLastLine('  x\n'), '  show(x)\n');
  assert.equal(autoShowLastLine('return x\n'), 'return x\n');
});

test('wrapNotebookExec and combinePythonOutput', () => {
  assert.equal(wrapNotebookExec('mymodel'), '_notebook_exec(' + JSON.stringify('show(mymodel)') + ')');
  assert.equal(combinePythonOutput(['line1', 'line2'], ['warn'], 'ret'), 'line1\nline2\nret\nwarn');
});

test('export helpers', () => {
  assert.equal(buildExportScript('stl'), [
    'if _last_shown_obj is None:',
    '    print("NO_MODEL_SHOWN_YET")',
    'else:',
    '    export(_last_shown_obj, "/tmp/download.stl")',
    '    print("EXPORT_OK")',
  ].join('\n'));
  assert.equal(exportPathForFormat('3mf'), '/tmp/download.3mf');
});

test('resolveWorkerResponse handles result and error', () => {
  {
    const pending = new Map([[
      7, {
        resolve: (v) => {
          resolved = v;
        },
        reject: () => {},
      }
    ]]);
    let resolved;
    assert.equal(resolveWorkerResponse(pending, {type: MSG.RESULT, id: 7, output: 'ok'}), true);
    assert.deepEqual(resolved, {type: MSG.RESULT, id: 7, output: 'ok'});
    assert.equal(pending.size, 0);
  } {const pending = new Map([[
       9, {
         resolve: () => {},
         reject: (e) => {
           rejected = e;
         },
       }
     ]]);
     let rejected; assert.equal(
       resolveWorkerResponse(pending, {type: MSG.ERROR, id: 9, message: 'boom', crashed: true}), true);
     assert.equal(rejected.message, 'boom'); assert.equal(rejected.crashed, true);}
});

class MockWorker
{
  constructor()
  {
    this.onmessage = null;
    this.postMessage = (msg) => {
      queueMicrotask(() => this._handle(msg));
    };
  }

  _reply(data)
  {
    if (typeof this.onmessage === 'function') this.onmessage({data});
  }

  _handle(msg)
  {
    if (msg.type === MSG.INIT) {
      this._reply({type: MSG.RESULT, id: msg.id, ready: true});
      return;
    }
    if (msg.type === MSG.EVALUATE) {
      this._reply({
        type: MSG.RESULT,
        id: msg.id,
        output: 'OUT:' + msg.code.slice(0, 8),
      });
      return;
    }
    if (msg.type === MSG.EXPORT) {
      this._reply({type: MSG.EXPORT_RESULT, id: msg.id, status: 'ok', format: msg.format});
    }
  }

  terminate() {}
}

test('NotebookKernel serializes requests against a mock worker', async () => {
  const kernel = new NotebookKernel({workerUrl: 'unused'});
  const worker = new MockWorker();
  kernel._attachWorker(worker);

  const initResult = await kernel._enqueue(() => kernel._request(MSG.INIT, {}));
  kernel.ready = Boolean(initResult.ready);
  assert.equal(kernel.ready, true);

  const evalResult = await kernel.evaluate('show(cube(1))');
  assert.equal(evalResult.output, 'OUT:show(cub');

  const exportResult = await kernel.exportModel('stl');
  assert.equal(exportResult.status, 'ok');

  let serial = '';
  const p1 = kernel.evaluate('first').then(r => {
    serial += '1:' + r.output;
  });
  const p2 = kernel.evaluate('second').then(r => {
    serial += '2:' + r.output;
  });
  await Promise.all([p1, p2]);
  assert.equal(serial, '1:OUT:first2:OUT:second');
});

test('NotebookKernel rejects in-flight and queued requests on termination', async () => {
  const kernel = new NotebookKernel({workerUrl: 'unused'});
  kernel._attachWorker({
    postMessage: () => {},
    terminate: () => {},
  });
  kernel.ready = true;

  const inFlight = kernel.evaluate('first');
  const queued = kernel.evaluate('second');
  kernel.terminate();

  await assert.rejects(inFlight, /Worker terminated/);
  await assert.rejects(queued, /Worker terminated/);
});
