import assert from 'node:assert/strict';

import {extractGeometryPayloads, mergeGeometries} from './notebook-geometry.mjs';

const blue = {
  type: '3d',
  vertices: [[0, 0, 0], [1, 0, 0], [0, 1, 0]],
  faces: [[0, 1, 2]],
  colors: [[0.2, 0.4, 0.6, 1]],
  color_indices: [0],
};
const yellow = {
  type: '3d',
  vertices: [[10, 0, 0], [11, 0, 0], [10, 1, 0]],
  faces: [[0, 1, 2]],
  colors: [[1, 0.9, 0.4, 1]],
  color_indices: [0],
};

const framed = [
  'before',
  '@@MIME_START@@' + JSON.stringify(blue) + '@@MIME_END@@',
  'between',
  '@@MIME_START@@' + JSON.stringify(yellow) + '@@MIME_END@@',
  'after',
].join('\n');
const extracted = extractGeometryPayloads(framed);
assert.deepEqual(extracted.geometries, [blue, yellow]);
assert.equal(extracted.text, 'before\n\nbetween\n\nafter');

const merged = mergeGeometries(extracted.geometries);
assert.deepEqual(merged.vertices, [...blue.vertices, ...yellow.vertices]);
assert.deepEqual(merged.faces, [[0, 1, 2], [3, 4, 5]]);
assert.deepEqual(merged.colors, [...blue.colors, ...yellow.colors]);
assert.deepEqual(merged.color_indices, [0, 1]);

const malformed = '@@MIME_START@@not-json@@MIME_END@@';
assert.deepEqual(extractGeometryPayloads(malformed), {
  geometries: [],
  text: malformed,
});
