import fs from 'node:fs';
import path from 'node:path';
import {pathToFileURL} from 'node:url';

function usage()
{
  console.error('Usage: node python-sandbox-runner.mjs --wasm-dir DIR --input FILE --output FILE');
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

if (!wasmDir || !inputFile || !outputFile) {
  usage();
  process.exit(2);
}

const modulePath = path.join(wasmDir, 'pythonscad.js');
const source = fs.readFileSync(inputFile, 'utf8');
const {default: OpenSCAD} = await import(pathToFileURL(modulePath).href);

const stdout = [];
const stderr = [];
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
    stderr.push(e?.stack || e?.message || String(e));
  }
}

if (exitCode !== 0) {
  for (const line of stdout) console.error(line);
  for (const line of stderr) console.error(line);
  process.exit(exitCode);
}

const csg = mod.FS.readFile('/work/output.csg', {encoding: 'utf8'});
fs.writeFileSync(outputFile, csg);
