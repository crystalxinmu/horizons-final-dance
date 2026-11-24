// 把源码与 three.js 打包成一个普通的 <script>，内联进 index.html，得到可离线双击打开的单文件。
import { build } from 'esbuild';
import { mkdir, readFile, writeFile } from 'node:fs/promises';

const result = await build({
  entryPoints: ['src/main.js'], bundle: true, minify: true, format: 'iife', target: ['es2020'],
  write: false, legalComments: 'none',
});
// 源码里若出现字面量 </script 会提前结束 HTML 的脚本标签
const js = result.outputFiles[0].text.replace(/<\/script/gi, '<\\/script');

const html = await readFile('index.html', 'utf8');
const importmap = /[ \t]*<script type="importmap">[\s\S]*?<\/script>\n/;
const entry = /<script type="module" src="\.\/src\/main\.js"><\/script>/;
if (!importmap.test(html) || !entry.test(html)) throw new Error('index.html 的结构不符合预期，无法内联');

const out = html.replace(importmap, '').replace(entry, () => `<script>\n${js}\n</script>`);
await mkdir('dist', { recursive: true });
await writeFile('dist/视界终舞.html', out);
await writeFile('dist/index.html', out);       // 内容相同；GitHub Pages 等托管服务以 index.html 为入口
console.log(`已生成 dist/视界终舞.html 与 dist/index.html（${(Buffer.byteLength(out) / 1024).toFixed(0)} KB，其中脚本 ${(Buffer.byteLength(js) / 1024).toFixed(0)} KB）`);
