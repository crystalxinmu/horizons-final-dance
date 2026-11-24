// 零依赖的静态服务器，用于开发版（index.html 通过 importmap 从 CDN 加载 three.js）
import http from 'node:http';
import { readFile, stat } from 'node:fs/promises';
import { extname, join, normalize, sep } from 'node:path';
import { fileURLToPath } from 'node:url';

const root = fileURLToPath(new URL('.', import.meta.url));
const port = Number(process.env.PORT ?? 5173);
const types = {
  '.html': 'text/html; charset=utf-8', '.js': 'text/javascript; charset=utf-8', '.css': 'text/css; charset=utf-8',
  '.json': 'application/json', '.md': 'text/markdown; charset=utf-8', '.png': 'image/png',
};

http.createServer(async (req, res) => {
  try {
    const pathname = decodeURIComponent(new URL(req.url, 'http://localhost').pathname);
    const file = normalize(join(root, pathname));
    if (!file.startsWith(root.endsWith(sep) ? root : root + sep)) throw Object.assign(new Error(), { code: 'FORBIDDEN' });
    let target = file;
    if ((await stat(target).catch(() => null))?.isDirectory()) target = join(target, 'index.html');
    const body = await readFile(target);
    res.writeHead(200, { 'content-type': types[extname(target)] ?? 'application/octet-stream', 'cache-control': 'no-store' });
    res.end(body);
  } catch {
    res.writeHead(404, { 'content-type': 'text/plain; charset=utf-8' });
    res.end('404 找不到文件');
  }
}).listen(port, () => console.log(`视界终舞 → http://localhost:${port}/`));
