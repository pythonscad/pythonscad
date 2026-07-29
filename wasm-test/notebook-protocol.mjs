/** Shared notebook kernel protocol helpers (main thread + worker). */

export const MSG = {
  INIT: 'init',
  EVALUATE: 'evaluate',
  EXPORT: 'export',
  RESULT: 'result',
  EXPORT_RESULT: 'export-result',
  ERROR: 'error',
};

export const REPL_PRELUDE = [
  'import ast as _ast',
  'import json as _json',
  'import pythonscad as _pyscad',
  'import openscad as _oscad',
  'from pythonscad import *',
  '_orig_show = _pyscad.show',
  '_last_shown_obj = None',
  'def _emit_geometry(obj):',
  '    global _last_shown_obj',
  '    _last_shown_obj = obj',
  '    try:',
  '        result = obj.mesh()',
  '        if len(result) > 0 and len(result[0]) > 0 and len(result[0][0]) == 2:',
  '            outlines_list = [[[float(p[0]), float(p[1])] for p in outline] for outline in result]',
  '            print(\'@@MIME_START@@\' + _json.dumps({\'type\': \'2d\', \'outlines\': outlines_list}) + \'@@MIME_END@@\')',
  '        else:',
  '            verts, faces, colors, color_indices = obj.mesh(color=True)',
  '            verts_list  = [[float(c) for c in v] for v in verts]',
  '            faces_list  = [[int(i) for i in f] for f in faces]',
  '            colors_list = [[float(x) for x in col] for col in colors]',
  '            cidx_list   = [int(i) for i in color_indices]',
  '            print(\'@@MIME_START@@\' + _json.dumps({',
  '                \'type\': \'3d\',',
  '                \'vertices\': verts_list,',
  '                \'faces\': faces_list,',
  '                \'colors\': colors_list,',
  '                \'color_indices\': cidx_list',
  '            }) + \'@@MIME_END@@\')',
  '    except Exception as _e:',
  '        print(\'GEOMETRY_ERROR: \' + repr(_e))',
  'def show(obj, *a, **kw):',
  '    r = _orig_show(obj, *a, **kw)',
  '    _emit_geometry(obj)',
  '    return r',
  'def _notebook_show(obj, *a, **kw):',
  '    r = obj.show(*a, **kw)',
  '    if isinstance(obj, _probe_base_class):',
  '        _emit_geometry(obj)',
  '    return r',
  'class _NotebookShowTransformer(_ast.NodeTransformer):',
  '    def visit_Call(self, node):',
  '        node = self.generic_visit(node)',
  '        if isinstance(node.func, _ast.Attribute) and node.func.attr == \'show\':',
  '            return _ast.copy_location(_ast.Call(',
  '                func=_ast.Name(id=\'_notebook_show\', ctx=_ast.Load()),',
  '                args=[node.func.value, *node.args],',
  '                keywords=node.keywords), node)',
  '        return node',
  'def _notebook_exec(source):',
  '    tree = _ast.parse(source, filename=\'<notebook>\', mode=\'exec\')',
  '    tree = _NotebookShowTransformer().visit(tree)',
  '    _ast.fix_missing_locations(tree)',
  '    exec(compile(tree, \'<notebook>\', \'exec\'), globals(), globals())',
  '_pyscad.show = show',
  '_oscad.show = show',
  '_probe_base_class = type(cube(1))',
  'print(\'PRELUDE_LOADED_OK\')',
].join('\n');

const RESERVED_LAST_LINE = new Set([
  'pass',
  'break',
  'continue',
  'return',
  'else',
  'elif',
  'try',
  'except',
  'finally',
]);

/** Jupyter-style: bare identifier on the last non-empty line becomes show(...). */
export function autoShowLastLine(code)
{
  const lines = String(code ?? '').split('\n');
  let lastIdx = lines.length - 1;
  while (lastIdx >= 0 && lines[lastIdx].trim() === '') lastIdx--;
  if (lastIdx < 0) return String(code ?? '');

  const originalLine = lines[lastIdx];
  const indent = originalLine.match(/^(\s*)/)[1];
  const lastLine = originalLine.trim();
  const isBareIdentifier = /^[A-Za-z_][A-Za-z0-9_]*$/.test(lastLine);
  if (isBareIdentifier && !RESERVED_LAST_LINE.has(lastLine)) {
    lines[lastIdx] = indent + 'show(' + lastLine + ')';
    return lines.join('\n');
  }
  return lines.join('\n');
}

export function wrapNotebookExec(code)
{
  const finalCode = autoShowLastLine(code);
  return '_notebook_exec(' + JSON.stringify(finalCode) + ')';
}

export function buildExportScript(format, path = '/tmp/download.' + format)
{
  const safePath = String(path);
  return [
    'if _last_shown_obj is None:',
    '    print("NO_MODEL_SHOWN_YET")',
    'else:',
    '    export(_last_shown_obj, ' + JSON.stringify(safePath) + ')',
    '    print("EXPORT_OK")',
  ].join('\n');
}

export function combinePythonOutput(stdoutLines, stderrLines, returnValue)
{
  const fromBuffers = (stdoutLines || []).join('\n');
  const fromStderr = (stderrLines || []).join('\n');
  return [fromBuffers, String(returnValue ?? '').trim(), fromStderr]
    .filter(s => s && s.trim())
    .join('\n');
}

export function exportPathForFormat(format)
{
  return '/tmp/download.' + String(format || 'stl');
}

/** Classify a worker->main response for a pending request id. */
export function resolveWorkerResponse(pending, data)
{
  const entry = pending.get(data.id);
  if (!entry) return false;

  pending.delete(data.id);
  if (data.type === MSG.ERROR) {
    const err = new Error(data.message || 'Worker error');
    err.crashed = Boolean(data.crashed);
    entry.reject(err);
    return true;
  }
  if (data.type === MSG.EXPORT_RESULT) {
    entry.resolve(data);
    return true;
  }
  if (data.type === MSG.RESULT) {
    entry.resolve(data);
    return true;
  }
  entry.reject(new Error('Unexpected worker message type: ' + data.type));
  return true;
}
