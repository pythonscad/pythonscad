import fs from 'node:fs';
import os from 'node:os';
import path from 'node:path';
import {pathToFileURL} from 'node:url';

function usage()
{
  console.error(
    'Usage: node python-sandbox-runner.mjs --wasm-dir DIR --input FILE --output FILE ' +
    '[--host-output-dir DIR --manifest FILE]');
}

function parseArgs(argv)
{
  const args = {};
  for (let i = 2; i < argv.length; i += 2) {
    const key = argv[i];
    const value = argv[i + 1];
    if (!key || !key.startsWith('--') || value === undefined) {
      usage();
      process.exit(2);
    }
    args[key.slice(2)] = value;
  }
  return args;
}

const args = parseArgs(process.argv);
const wasmDir = args['wasm-dir'];
const inputFile = args.input;
const outputFile = args.output;
const hostOutputDir = args['host-output-dir'];
const manifestFile = args.manifest;

if (!wasmDir || !inputFile || !outputFile) {
  usage();
  process.exit(2);
}

const modulePath = path.join(wasmDir, 'pythonscad.js');
const source = fs.readFileSync(inputFile, 'utf8');

function removeTreeSync(dir)
{
  if (typeof fs.rmSync === 'function') {
    fs.rmSync(dir, {recursive: true, force: true});
  } else {
    fs.rmdirSync(dir, {recursive: true});
  }
}

async function importOpenSCADModule()
{
  try {
    return (await import(pathToFileURL(modulePath).href)).default;
  } catch (error) {
    if (!(error instanceof SyntaxError)) throw error;
    const tempDir = fs.mkdtempSync(path.join(os.tmpdir(), 'pythonscad-wasm-module-'));
    const tempModulePath = path.join(tempDir, 'pythonscad.mjs');
    try {
      fs.copyFileSync(modulePath, tempModulePath);
      return (await import(pathToFileURL(tempModulePath).href)).default;
    } finally {
      removeTreeSync(tempDir);
    }
  }
}

const RESERVED_WINDOWS_NAMES = new Set([
  'con',  'prn',  'aux',  'nul',  'com1', 'com2', 'com3', 'com4', 'com5', 'com6', 'com7',
  'com8', 'com9', 'lpt1', 'lpt2', 'lpt3', 'lpt4', 'lpt5', 'lpt6', 'lpt7', 'lpt8', 'lpt9',
]);

function safeRelativePath(relativePath)
{
  if (!relativePath || relativePath.includes('\0') || relativePath.startsWith('/')) return null;
  if (/[\t\r\n]/.test(relativePath)) return null;
  if (relativePath.includes(':')) return null;
  if (/^[A-Za-z]:/.test(relativePath) || relativePath.startsWith('\\\\')) return null;

  const normalized = path.posix.normalize(relativePath.replace(/\\/g, '/'));
  if (!normalized || normalized === '.' || normalized.startsWith('../') || normalized === '..')
    return null;

  const parts = normalized.split('/');
  for (const part of parts) {
    if (!part || part === '.' || part === '..') return null;
    const stem = part.split('.')[0].toLowerCase();
    if (RESERVED_WINDOWS_NAMES.has(stem)) return null;
  }
  return normalized;
}

function listFilesRecursive(fsApi, root)
{
  try {
    const stat = fsApi.stat(root);
    if (!fsApi.isDir(stat.mode)) return [];
  } catch {
    return [];
  }

  const files = [];
  const walk = (dir, relPrefix) => {
    for (const name of fsApi.readdir(dir)) {
      if (name === '.' || name === '..') continue;
      const virtualPath = path.posix.join(dir, name);
      const relPath = relPrefix ? path.posix.join(relPrefix, name) : name;
      const stat = fsApi.stat(virtualPath);
      if (fsApi.isDir(stat.mode)) walk(virtualPath, relPath);
      else if (fsApi.isFile(stat.mode)) files.push({virtualPath, relPath});
    }
  };
  walk(root, '');
  return files;
}

function collectSandboxOutputs(fsApi)
{
  if (!hostOutputDir || !manifestFile) return;

  const roots = ['/work/out', '/outputs'];
  const copied = new Set();
  const manifest = [];

  for (const root of roots) {
    for (const file of listFilesRecursive(fsApi, root)) {
      const relPath = safeRelativePath(file.relPath);
      if (!relPath || copied.has(relPath)) continue;
      copied.add(relPath);

      const data = fsApi.readFile(file.virtualPath);
      const hostPath = path.join(hostOutputDir, ...relPath.split('/'));
      fs.mkdirSync(path.dirname(hostPath), {recursive: true});
      fs.writeFileSync(hostPath, Buffer.from(data));
      manifest.push(`${relPath}\t${hostPath}\t${data.length}`);
    }
  }

  fs.writeFileSync(manifestFile, manifest.join('\n') + (manifest.length ? '\n' : ''));
}

const stdout = [];
const stderr = [];

async function main()
{
  const OpenSCAD = await importOpenSCADModule();
  const mod = await OpenSCAD({
    locateFile: (name) => path.join(wasmDir, name),
    print: (line) => stdout.push(line),
    printErr: (line) => stderr.push(line),
  });

  mod.FS.mkdir('/work');
  mod.FS.writeFile('/work/input.py', source);

  let exitCode = 0;
  try {
    mod.callMain(['-o', '/work/output.csg', '--python=native', '/work/input.py']);
  } catch (e) {
    if (e && typeof e.status === 'number') {
      exitCode = e.status;
    } else {
      exitCode = 1;
      stderr.push((e && (e.stack || e.message)) || String(e));
    }
  }

  if (exitCode !== 0) {
    for (const line of stdout) console.error(line);
    for (const line of stderr) console.error(line);
    process.exit(exitCode);
  }

  const csg = mod.FS.readFile('/work/output.csg', {encoding: 'utf8'});
  fs.writeFileSync(outputFile, csg);
  collectSandboxOutputs(mod.FS);
}

main().catch((e) => {
  console.error((e && (e.stack || e.message)) || String(e));
  process.exit(1);
});
