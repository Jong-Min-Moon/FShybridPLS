/**
 * Capture animation/index.html into hybridpls-demo.gif for the GitHub README.
 * Usage: node record_gif.mjs
 */
import { createServer } from "node:http";
import { readFile } from "node:fs/promises";
import { createWriteStream } from "node:fs";
import path from "node:path";
import { fileURLToPath } from "node:url";
import { chromium } from "playwright";
import { spawnSync } from "node:child_process";

const __dirname = path.dirname(fileURLToPath(import.meta.url));
const ROOT = __dirname;
const OUT_GIF = path.join(ROOT, "hybridpls-demo.gif");
const FRAME_DIR = path.join(ROOT, "_frames");

const MIME = {
  ".html": "text/html; charset=utf-8",
  ".js": "text/javascript; charset=utf-8",
  ".json": "application/json",
  ".css": "text/css",
  ".png": "image/png",
  ".svg": "image/svg+xml"
};

function startServer() {
  const server = createServer(async (req, res) => {
    try {
      const url = new URL(req.url || "/", "http://127.0.0.1");
      let rel = decodeURIComponent(url.pathname);
      if (rel === "/") rel = "/index.html";
      const file = path.join(ROOT, rel.replace(/^\/+/, ""));
      if (!file.startsWith(ROOT)) {
        res.writeHead(403); res.end("forbidden"); return;
      }
      const data = await readFile(file);
      res.writeHead(200, { "Content-Type": MIME[path.extname(file)] || "application/octet-stream" });
      res.end(data);
    } catch {
      res.writeHead(404); res.end("not found");
    }
  });
  return new Promise((resolve) => {
    server.listen(0, "127.0.0.1", () => resolve(server));
  });
}

async function main() {
  const { mkdirSync, rmSync, existsSync } = await import("node:fs");
  if (existsSync(FRAME_DIR)) rmSync(FRAME_DIR, { recursive: true });
  mkdirSync(FRAME_DIR, { recursive: true });

  const server = await startServer();
  const { port } = server.address();
  const url = `http://127.0.0.1:${port}/index.html?record=1`;

  console.log("Serving", url);
  const browser = await chromium.launch({ headless: true });
  const page = await browser.newPage({
    viewport: { width: 860, height: 560 },
    deviceScaleFactor: 1
  });
  await page.goto(url, { waitUntil: "networkidle" });
  await page.waitForFunction(() => window.__ANIM__ && window.__ANIM__.ready);

  const totalMs = await page.evaluate(() => window.__ANIM__.totalMs);
  // ~5 fps over the full timeline → ~15s GIF playback of the whole story
  const stepMs = 200;
  const frames = [];
  console.log(`Capturing ${Math.floor(totalMs / stepMs) + 1} frames over ${totalMs}ms…`);

  for (let t = 0, i = 0; t <= totalMs; t += stepMs, i++) {
    await page.evaluate((ms) => window.__ANIM__.seek(ms), Math.min(t, totalMs));
    // Let SVG paint; freeze time-based pulses at a stable phase by seeking again
    await page.waitForTimeout(30);
    const file = path.join(FRAME_DIR, `f${String(i).padStart(4, "0")}.png`);
    await page.locator(".frame").screenshot({ path: file, type: "png" });
    frames.push(file);
    if (i % 25 === 0) console.log(`  frame ${i} @ ${t}ms`);
  }

  await browser.close();
  server.close();

  console.log("Assembling GIF with Python/Pillow…");
  const py = `
from PIL import Image
import os
frames = ${JSON.stringify(frames).replace(/\\/g, "/")}
out = r"""${OUT_GIF.replace(/\\/g, "/")}"""
imgs = []
for f in frames:
    im = Image.open(f).convert("P", palette=Image.ADAPTIVE, colors=128)
    imgs.append(im)
imgs[0].save(
    out,
    save_all=True,
    append_images=imgs[1:],
    duration=120,
    loop=0,
    optimize=True,
    disposal=2
)
print("wrote", out, "frames", len(imgs), "bytes", os.path.getsize(out))
`;
  const r = spawnSync("python", ["-c", py], { encoding: "utf-8" });
  if (r.status !== 0) {
    console.error(r.stdout, r.stderr);
    process.exit(1);
  }
  console.log(r.stdout.trim());
  rmSync(FRAME_DIR, { recursive: true });
  console.log("Done:", OUT_GIF);
}

main().catch((e) => {
  console.error(e);
  process.exit(1);
});
