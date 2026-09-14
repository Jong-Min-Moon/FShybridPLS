# Hybrid PLS animation

Educational walkthrough of Hybrid PLS (motivation → associations → hybrid \(W\) → PLS weights → one score \(\rho\)).

## View

- Open [`index.html`](index.html) in a browser (loads sibling `animation_data.js`).
- README on GitHub shows the looping demo GIF [`hybridpls-demo.gif`](hybridpls-demo.gif) (GitHub cannot run interactive HTML inside README.md).

## Live demo (GitHub Pages)

The Actions-based Pages workflow was removed (it fails until Pages exists, and `GITHUB_TOKEN` cannot create Pages).

Use **branch deploy** instead:

1. Commit/push the `docs/` folder at the repo root (copy of this animation site).
2. Open https://github.com/Jong-Min-Moon/FShybridPLS/settings/pages
3. **Build and deployment → Source:** Deploy from a branch
4. Branch: `main` · folder: `/docs` · Save

Site: https://jong-min-moon.github.io/FShybridPLS/

Until then, use the CDN mirror in the main README.

## Record / refresh the GIF

```bash
cd animation
npm install
npx playwright install chromium
npm run record-gif
```

This writes `hybridpls-demo.gif`. After regenerating, also copy it into `docs/` if you use Pages:

```bash
copy hybridpls-demo.gif ..\docs\
copy index.html ..\docs\
copy animation_data.js ..\docs\
```
