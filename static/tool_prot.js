// tool_prot.js

// === Data placeholders ===
// protData is injected by Python as: {gene: {plddt:[], fpocket:[], pdb:"..."}}
var currentGene = null;
let viewer = null;

// === Viridis Color Map ===
function viridisColor(value) {
  const viridis = [
    [68, 1, 84], [71, 44, 122], [59, 81, 139], [44, 113, 142],
    [33, 144, 141], [39, 173, 129], [92, 200, 99], [170, 220, 50],
    [253, 231, 37]
  ];
  value = Math.min(1, Math.max(0, value));
  const idx = value * (viridis.length - 1);
  const i = Math.floor(idx);
  const f = idx - i;

  if (i >= viridis.length - 1) {
    const [r, g, b] = viridis[viridis.length - 1];
    return (r << 16) | (g << 8) | b;
  }

  const rgb = [0, 1, 2].map(j =>
    Math.round(viridis[i][j] + f * (viridis[i + 1][j] - viridis[i][j]))
  );
  return (rgb[0] << 16) | (rgb[1] << 8) | rgb[2];
}

// === Legend Drawing ===
function drawLegend() {
  const canvas = document.getElementById("colorLegend");
  if (!canvas) return;
  const ctx = canvas.getContext("2d");
  const width = canvas.width, height = canvas.height;
  for (let i = 0; i < width; i++) {
    const t = i / (width - 1);
    const rgbInt = viridisColor(t);
    const hex = rgbInt.toString(16).padStart(6, "0");
    ctx.fillStyle = "#" + hex;
    ctx.fillRect(i, 0, 1, height);
  }
}

// === Contribution Normalizers ===
function getPLDDTContributions(gene) {
  if (!protData[gene]) return [];
  return protData[gene].plddt || [];
}

function getFpocketContributions(gene) {
  if (!protData[gene]) return [];
  return protData[gene].fpocket || [];
}

// === Surface Coloring Helper ===
function colorSurface(contributions) {
  const scoreMap = {};
  contributions.forEach(c => {
    if (c.residue_no != null) scoreMap[c.residue_no] = c.score;
  });
  viewer.addSurface($3Dmol.SurfaceType.VDW, {
    opacity: 1.0,
    colorfunc: atom => {
      const s = scoreMap[atom.resi];
      return s !== undefined ? viridisColor(s) : 0xAAAAAA;
    }
  });
  document.getElementById("legend").style.display = "flex";
  drawLegend();
}

// === Apply Selected Style ===
function applyStyle(gene, style) {
  if (!protData[gene]) return;

  viewer.clear();
  viewer.addModel(protData[gene].pdb, "pdb");

  if (style === "cartoon") {
    viewer.setStyle({}, { cartoon: { color: "spectrum" } });
    document.getElementById("legend").style.display = "none";
  }
  else if (style === "surface-plddt") {
    colorSurface(getPLDDTContributions(gene));
  }
  else if (style === "surface-fpocket") {
    colorSurface(getFpocketContributions(gene));
  }

  viewer.zoomTo();
  viewer.render();
}

// === Init on Page Load ===
document.addEventListener("DOMContentLoaded", () => {
  const cssBg = getComputedStyle(document.body).getPropertyValue("--bg-panel").trim() || "#ffffff";
  viewer = $3Dmol.createViewer("viewer-container", { backgroundColor: cssBg });

  const protSelect = document.getElementById("proteinSelect");

  // Populate protein dropdown from injected protData
  Object.keys(protData).forEach(gene => {
    const opt = document.createElement("option");
    opt.value = gene;
    opt.textContent = gene;
    protSelect.appendChild(opt);
  });

  // Default: first gene if available
  if (Object.keys(protData).length > 0) {
    currentGene = Object.keys(protData)[0];
    protSelect.value = currentGene;
    applyStyle(currentGene, "cartoon");
  }

  // Style change handler
  document.getElementById("styleSelect").addEventListener("change", e => {
    if (currentGene) applyStyle(currentGene, e.target.value);
  });

  // Protein change handler
  protSelect.addEventListener("change", e => {
    currentGene = e.target.value;
    const style = document.getElementById("styleSelect").value;
    applyStyle(currentGene, style);
  });

  // === Theme sync via postMessage ===
  window.addEventListener("message", (event) => {
  if (!event.data || event.data.type !== "theme") return;

  const theme = event.data.theme;

  // Toggle the body class so CSS variables update
  if (theme === "light") {
    document.body.classList.add("light");
  } else {
    document.body.classList.remove("light");
  }

  // Update 3Dmol viewer background
  if (viewer) {
    const bg = theme === "light" ? "#ffffff" : "#0e1430";
    viewer.setBackgroundColor(bg);
    viewer.render();
  }
 });

});


