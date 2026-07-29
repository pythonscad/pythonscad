const MARK_START = '@@MIME_START@@';
const MARK_END = '@@MIME_END@@';

export function extractGeometryPayloads(output)
{
  const text = String(output ?? '');
  const geometries = [];
  const textParts = [];
  let cursor = 0;

  while (cursor < text.length) {
    const start = text.indexOf(MARK_START, cursor);
    if (start === -1) {
      textParts.push(text.slice(cursor));
      break;
    }
    textParts.push(text.slice(cursor, start));
    const end = text.indexOf(MARK_END, start + MARK_START.length);
    if (end === -1) {
      textParts.push(text.slice(start));
      break;
    }
    const jsonStr = text.slice(start + MARK_START.length, end);
    try {
      geometries.push(JSON.parse(jsonStr));
    } catch (e) {
      textParts.push(text.slice(start, end + MARK_END.length));
    }
    cursor = end + MARK_END.length;
  }

  return {geometries, text: textParts.join('').trim()};
}

export function mergeGeometries(geometries)
{
  const merged = {type: '3d', vertices: [], faces: [], colors: [], color_indices: []};

  for (const geometry of geometries) {
    if (!geometry) continue;
    const vertices = Array.isArray(geometry.vertices) ? geometry.vertices : [];
    const faces = Array.isArray(geometry.faces) ? geometry.faces : [];
    const colors = Array.isArray(geometry.colors) ? geometry.colors : [];
    const colorIndices = Array.isArray(geometry.color_indices) ? geometry.color_indices : [];
    const vertexOffset = merged.vertices.length;
    const colorOffset = merged.colors.length;

    merged.vertices.push(...vertices);
    merged.colors.push(...colors);
    for (let i = 0; i < faces.length; ++i) {
      merged.faces.push(faces[i].map(index => index + vertexOffset));
      const colorIndex = i < colorIndices.length ? colorIndices[i] : -1;
      merged.color_indices.push(
        colorIndex >= 0 && colorIndex < colors.length ? colorIndex + colorOffset : -1);
    }
  }

  return merged;
}
