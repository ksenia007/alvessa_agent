// Variant summary card
// Mirrors structure used by static/gene_summary.js

const MAX_VARIANT_ROWS = 400;

export async function renderVariantSummary(st, deps = {}) {
  const byId      = deps.byId || ((id) => document.getElementById(id));
  const show      = deps.show || ((id) => { const el = byId(id); if (el) el.style.display = el.classList.contains("grid") ? "grid" : "block"; });
  const hide      = deps.hide || ((id) => { const el = byId(id); if (el) el.style.display = "none"; });
  const fetchText = deps.fetchTextCached || (async (p) => (await fetch(p, { cache: "no-store" })).text());
  const parseCSV  = deps.parseCSV || defaultParseCSV;
  const esc       = deps.escapeHtml || defaultEscapeHtml;

  const mount = byId("variantSummary");
  const card  = byId("variantSummaryCard");
  if (!mount || !card) return;

  const variantsCsv = await loadCsvSafe("variants/variants.csv");
  if (!variantsCsv.length) {
    mount.innerHTML = "";
    hide("variantSummaryCard");
    return;
  }

  const baseOrder = [];
  const variantMeta = new Map();
  for (const row of variantsCsv) {
    const vid = normVariant(row.variant_id);
    if (!vid) continue;
    if (!baseOrder.includes(vid)) baseOrder.push(vid);
    variantMeta.set(vid, {
      variant: vid,
      chrom: row.primary_chrom || row.chrom || "",
      pos: row.primary_pos || row.pos || "",
      build: row.primary_build || row.build || "",
      organism: row.organism || "",
    });
  }

  const expectoStore = new Map();
  const alleleStore = new Map();
  const traitStore = new Map();

  const { entries, toolNames, hasContext } = await buildVariantGeneEntries({
    fetchText,
    parseCSV,
  });
  const alleleData = await loadAlleleFrequencies({ fetchText, parseCSV });
  const traitData  = await loadTraitData({ fetchText, parseCSV });

  const variantIds = baseOrder.length ? baseOrder : uniquePreserve([...entries.keys()]);
  const rows = [];
  for (const vid of variantIds) {
    const genes = entries.get(vid);
    if (!genes || !genes.size) continue;
    const meta = variantMeta.get(vid) || {};
    const chrom = meta.chrom || meta.primary_chrom || "";
    const pos = meta.pos || meta.primary_pos || "";
    const locationLabel = formatLocationLabel(chrom, pos);
    const geneList = Array.from(genes.values()).sort((a, b) => a.label.localeCompare(b.label));
    for (const entry of geneList) {
      const cells = [
        vid,
        locationLabel,
      ];
      cells.push(renderAlleleCell(vid, locationLabel, alleleData.get(vid), alleleStore));
      cells.push(renderTraitCell(vid, entry.geneKey, traitData, traitStore, entry.label));
      if (hasContext) cells.push(displayValue(entry.context));
      for (const tool of toolNames) {
        cells.push(formatPredictionCell(
          tool,
          entry.predictions.get(tool),
          expectoStore,
          {
            variantId: vid,
            variantLabel: vid,
            locationLabel,
            geneLabel: entry.label || entry.geneKey || "—",
          }
        ));
      }
      rows.push({
        gene: entry.label || entry.geneKey || "—",
        cells,
      });
    }
  }

  if (!rows.length) {
    mount.innerHTML = `<div class="muted">No variant-level gene annotations available.</div>`;
    show("variantSummaryCard");
    return;
  }

  const header = [];
  const columnStyles = [];
  header.push("Variant");
  columnStyles.push("min-width:160px; white-space:nowrap; font-family:var(--mono, ui-monospace);");
  header.push("Location");
  columnStyles.push("min-width:190px; white-space:nowrap;");
  header.push("TopMed AF");
  columnStyles.push("min-width:120px; white-space:nowrap;");

  header.push("Top Traits");
  columnStyles.push("min-width:200px; white-space:normal;");

  if (hasContext) {
    header.push("Context");
    columnStyles.push("min-width:160px; white-space:nowrap;");
  }
  for (const tool of toolNames) {
    header.push(toolLabel(tool));
    columnStyles.push("min-width:150px;");
  }

  const limitedRows = rows.slice(0, MAX_VARIANT_ROWS);
  const byGene = new Map();
  for (const row of limitedRows) {
    if (!byGene.has(row.gene)) byGene.set(row.gene, []);
    byGene.get(row.gene).push(row.cells);
  }
  const geneSections = Array.from(byGene.entries()).sort((a, b) => a[0].localeCompare(b[0]));
  const tablesHtml = geneSections.map(([gene, geneRows]) => `
    <div style="margin-top:12px;">
      <div class="muted" style="font-weight:600; border-top:1px dashed var(--line); padding-top:8px;">
        ${esc(gene)} (${geneRows.length})
      </div>
      <div style="margin-top:6px;">
        ${renderMiniTable(header, geneRows, columnStyles)}
      </div>
    </div>
  `).join("");

  const extraNote = rows.length > MAX_VARIANT_ROWS
    ? `<div class="muted" style="margin-top:6px;">Showing first ${MAX_VARIANT_ROWS} of ${rows.length} variant-gene rows.</div>`
    : "";
  const toolsList = toolNames.length ? toolNames.map(toolLabel).join(", ") : "—";

  mount.innerHTML = `
    <div class="muted" style="margin-bottom:8px;">
      Variant-gene rows: ${rows.length}. Genes rendered: ${geneSections.length}. Tools detected: ${esc(toolsList)}.
    </div>
    ${tablesHtml || `<div class="muted">No gene-level variant rows to display.</div>`}
    ${extraNote}
  `;
  show("variantSummaryCard");
  wireExpectoPanels(mount, expectoStore);
  wireAllelePanels(mount, alleleStore);
  wireTraitPanels(mount, traitStore);

  async function loadCsvSafe(path) {
    try {
      const txt = await fetchText(path);
      return parseCSV(txt, ",");
    } catch {
      return [];
    }
  }
}

async function buildVariantGeneEntries({ fetchText, parseCSV }) {
  const entries = new Map();
  const toolSet = new Set();
  let hasContext = false;

  const contextRows = await loadCsv("variants/per_gene_context.csv");
  for (const row of contextRows) {
    const vid = normVariant(row.variant_id);
    const geneKey = normGene(row.gene);
    if (!vid || !geneKey) continue;
    const entry = ensureEntry(entries, vid, geneKey, row.gene);
    if (row.context) {
      entry.context = row.context;
      hasContext = true;
    }
  }

  const predictionRows = await loadCsv("variants/functional_predictions.csv");
  for (const row of predictionRows) {
    const vid = normVariant(row.variant_id);
    const geneKey = normGene(row.gene);
    const tool = String(row.tool || "").trim();
    if (!vid || !geneKey || !tool) continue;
    const entry = ensureEntry(entries, vid, geneKey, row.gene);
    const score = (row.score || "").trim();
    entry.predictions.set(tool, score || "—");
    toolSet.add(tool);
  }

  return {
    entries,
    toolNames: Array.from(toolSet).sort((a, b) => a.localeCompare(b)),
    hasContext,
  };

  async function loadCsv(path) {
    try {
      const txt = await fetchText(path);
      return parseCSV(txt, ",");
    } catch {
      return [];
    }
  }
}

function ensureEntry(entries, variantId, geneKey, label) {
  if (!entries.has(variantId)) entries.set(variantId, new Map());
  const genes = entries.get(variantId);
  const key = geneKey || label || "—";
  if (!genes.has(key)) {
    genes.set(key, {
      geneKey: key,
      label: label ? label.trim() : key,
      context: "",
      predictions: new Map(),
    });
  }
  const entry = genes.get(key);
  if (label) entry.label = label.trim() || entry.label;
  return entry;
}

async function loadAlleleFrequencies({ fetchText, parseCSV }) {
  const rows = await (async () => {
    try {
      const txt = await fetchText("variants/allele_frequencies.csv");
      return parseCSV(txt, ",");
    } catch {
      return [];
    }
  })();
  const map = new Map();
  for (const row of rows) {
    const vid = normVariant(row.variant_id);
    if (!vid) continue;
    if (!map.has(vid)) map.set(vid, { rows: [], topmed: null });
    const bucket = map.get(vid);
    const raw = safeJson(row.raw_json) || {};
    const study = raw.study_name || row.source || row.population || "Study";
    const alleleFreq = parseNumber(row.af) ?? parseNumber(raw.allele_frequency);
    const alleleCount = parseNumber(row.ac) ?? parseNumber(raw.allele_count);
    const totalCount  = parseNumber(row.an) ?? parseNumber(raw.total_count);
    const obs = (() => {
      const rawObs = raw.observations ?? raw.observation;
      if (Array.isArray(rawObs)) return rawObs[0] || {};
      return rawObs || {};
    })();
    const ref = obs.deleted_sequence || row.ref || "";
    const alt = obs.inserted_sequence || row.alt || "";
    const keepRow = !ref || !alt || ref !== alt;
    if (keepRow) {
      bucket.rows.push({
        study,
        allele_frequency: alleleFreq,
        allele_count: alleleCount,
        total_count: totalCount,
        population: row.population || raw.population || "",
        source: row.source || raw.source || "",
        ref,
        alt,
      });
      if (bucket.topmed == null && String(study || "").toUpperCase() === "TOPMED" && Number.isFinite(alleleFreq)) {
        bucket.topmed = alleleFreq;
      }
    }
  }
  return map;
}

async function loadTraitData({ fetchText, parseCSV }) {
  const rows = await (async () => {
    try {
      const txt = await fetchText("variants/per_gene_traits.csv");
      return parseCSV(txt, ",");
    } catch {
      return [];
    }
  })();
  const map = new Map();
  for (const row of rows) {
    const vid = normVariant(row.variant_id);
    const geneKey = normGene(row.gene);
    if (!vid || !geneKey) continue;
    if (!map.has(vid)) map.set(vid, new Map());
    const geneMap = map.get(vid);
    if (!geneMap.has(geneKey)) geneMap.set(geneKey, []);
    const trait = row.trait || row.trait_name || "";
    const pVal = parseNumber(row.p_value);
    const risk = parseNumber(row.risk_score);
    geneMap.get(geneKey).push({
      trait,
      p_value: pVal,
      risk_score: risk,
    });
  }
  return map;
}

function renderAlleleCell(variantId, locationLabel, data, store) {
  if (!data || !Array.isArray(data.rows) || data.rows.length === 0) {
    return "—";
  }
  const filtered = dedupeAlleleRows(
    data.rows.filter((rec) => !rec.ref || !rec.alt || rec.ref !== rec.alt)
  );
  if (!filtered.length) return "—";

  const topmed = Number.isFinite(data.topmed) ? formatPercent(data.topmed) : "—";
  const id = `af_${rand7()}`;
  store.set(id, {
    variant: variantId,
    location: locationLabel,
    rows: filtered,
  });
  return {
    __html: `
      <div style="display:flex;flex-direction:column;gap:4px;">
        <div>${escapeInline(topmed)}</div>
        <button type="button" class="vs-link" data-af-id="${escapeInline(id)}">All studies</button>
      </div>
    `,
  };
}

function renderTraitCell(variantId, geneKey, traitData, store, geneLabel = "") {
  const geneMap = traitData.get(variantId);
  const traits = geneMap?.get(normGene(geneKey)) || [];
  if (!traits.length) return "—";
  const sorted = sortTraits(traits);
  const topTwo = sorted.slice(0, 2);
  const summary = topTwo
    .map((t) => `${escapeInline(t.trait || "Trait")} (${formatRisk(t.risk_score)}, p=${formatPValue(t.p_value)})`)
    .join("<br>");
  const id = `trait_${rand7()}`;
  store.set(id, {
    variant: variantId,
    gene: geneLabel || geneKey,
    rows: sorted,
  });
  return {
    __html: `
      <div style="display:flex;flex-direction:column;gap:4px;">
        <div>${summary || "—"}</div>
        <button type="button" class="vs-link" data-trait-id="${escapeInline(id)}">All traits</button>
      </div>
    `,
  };
}

function wireAllelePanels(root, store) {
  if (!store || !store.size) return;
  root.querySelectorAll("[data-af-id]").forEach((btn) => {
    const id = btn.getAttribute("data-af-id");
    const payload = store.get(id);
    if (!id || !payload) return;
    btn.addEventListener("click", (e) => {
      e.preventDefault();
      openAlleleModal(payload);
    });
  });
}

function wireTraitPanels(root, store) {
  if (!store || !store.size) return;
  root.querySelectorAll("[data-trait-id]").forEach((btn) => {
    const id = btn.getAttribute("data-trait-id");
    const payload = store.get(id);
    if (!id || !payload) return;
    btn.addEventListener("click", (e) => {
      e.preventDefault();
      openTraitModal(payload);
    });
  });
}

function openAlleleModal(data) {
  const subtitle = data.location ? `Location: ${data.location}` : "";
  const { body } = showVariantModal({
    title: `${data.variant} — allele frequencies`,
    subtitle,
    width: 720,
  });

  const controls = document.createElement("div");
  controls.style.display = "flex";
  controls.style.gap = "8px";
  controls.style.flexWrap = "wrap";

  const sortButtons = [
    { key: "total", label: "Sort by total count" },
    { key: "af", label: "Sort by allele frequency" },
  ].map(({ key, label }) => {
    const btn = document.createElement("button");
    btn.type = "button";
    btn.textContent = label;
    btn.className = "vs-link";
    btn.dataset.sortKey = key;
    controls.appendChild(btn);
    return btn;
  });

  const tableWrap = document.createElement("div");
  tableWrap.style.marginTop = "10px";

  body.appendChild(controls);
  body.appendChild(tableWrap);

  let sortKey = "total";
  const render = () => {
    tableWrap.innerHTML = renderAlleleTable(data.rows, sortKey);
    sortButtons.forEach((btn) => btn.classList.toggle("active", btn.dataset.sortKey === sortKey));
  };

  sortButtons.forEach((btn) => {
    btn.addEventListener("click", () => {
      sortKey = btn.dataset.sortKey || "total";
      render();
    });
  });

  render();
}

function renderAlleleTable(rows, sortKey) {
  const sorted = [...rows].sort((a, b) => {
    if (sortKey === "af") return (b.allele_frequency ?? -Infinity) - (a.allele_frequency ?? -Infinity);
    return (b.total_count ?? -Infinity) - (a.total_count ?? -Infinity);
  });
  const table = sorted
    .map((rec) => `
      <tr>
        <td class="wrap">${escapeInline(rec.study || "Study")}</td>
        <td class="wrap" style="text-align:center;">${escapeInline(rec.ref || "—")}</td>
        <td class="wrap" style="text-align:center;">${escapeInline(rec.alt || "—")}</td>
        <td class="wrap" style="text-align:right;">${formatPercent(rec.allele_frequency)}</td>
        <td class="wrap" style="text-align:right;">${formatInteger(rec.total_count)}</td>
      </tr>
    `)
    .join("");
  return `
    <div style="max-height:360px; overflow:auto;">
      <table class="vs-table">
        <thead>
          <tr>
            <th class="wrap">Study</th>
            <th class="wrap" style="text-align:center;">Deleted</th>
            <th class="wrap" style="text-align:center;">Inserted</th>
            <th class="wrap" style="text-align:right;">Allele Frequency</th>
            <th class="wrap" style="text-align:right;">Total Count</th>
          </tr>
        </thead>
        <tbody>${table}</tbody>
      </table>
    </div>
  `;
}

function openTraitModal(data) {
  const { body } = showVariantModal({
    title: `${data.variant} — ${data.gene} traits`,
    subtitle: "",
    width: 720,
  });
  const table = data.rows
    .map((rec) => `
      <tr>
        <td class="wrap">${escapeInline(rec.trait || "Trait")}</td>
        <td class="wrap" style="text-align:right;">${formatRisk(rec.risk_score)}</td>
        <td class="wrap" style="text-align:right;">${formatPValue(rec.p_value)}</td>
      </tr>
    `)
    .join("");
  body.innerHTML = `
    <div style="max-height:420px; overflow:auto;">
      <table class="vs-table">
        <thead>
          <tr>
            <th class="wrap">Trait</th>
            <th class="wrap" style="text-align:right;">Risk Score</th>
            <th class="wrap" style="text-align:right;">p-value</th>
          </tr>
        </thead>
        <tbody>${table}</tbody>
      </table>
    </div>
  `;
}

function showVariantModal({ title, subtitle, width = 680 }) {
  const overlay = document.createElement("div");
  overlay.className = "variant-overlay";
  overlay.style.position = "fixed";
  overlay.style.inset = "0";
  overlay.style.background = "rgba(2,4,12,0.65)";
  overlay.style.backdropFilter = "blur(4px)";
  overlay.style.zIndex = "9999";
  overlay.style.display = "flex";
  overlay.style.alignItems = "center";
  overlay.style.justifyContent = "center";

  const modal = document.createElement("div");
  modal.className = "variant-modal";
  modal.style.width = "90%";
  modal.style.maxWidth = `${width}px`;
  modal.style.maxHeight = "85vh";
  modal.style.overflow = "auto";
  modal.style.background = "var(--card)";
  modal.style.border = "1px solid var(--line)";
  modal.style.borderRadius = "18px";
  modal.style.boxShadow = "0 30px 80px rgba(0,0,0,0.45)";
  modal.style.padding = "18px 20px";

  const header = document.createElement("div");
  header.style.display = "flex";
  header.style.alignItems = "center";
  header.style.justifyContent = "space-between";
  header.style.gap = "16px";
  header.innerHTML = `
    <div>
      <div style="font-weight:700;">${escapeInline(title || "")}</div>
      ${subtitle ? `<div class="muted" style="font-size:13px;">${escapeInline(subtitle)}</div>` : ""}
    </div>
    <button type="button" class="variant-close" style="border:none;background:transparent;font-size:20px;cursor:pointer;color:var(--fg);">×</button>
  `;

  const body = document.createElement("div");
  body.className = "variant-modal-body";
  body.style.marginTop = "14px";

  modal.appendChild(header);
  modal.appendChild(body);
  overlay.appendChild(modal);
  document.body.appendChild(overlay);

  const cleanup = () => {
    overlay.remove();
    document.removeEventListener("keydown", escHandler);
  };
  const escHandler = (e) => {
    if (e.key === "Escape") cleanup();
  };

  overlay.addEventListener("click", (e) => {
    if (e.target === overlay) cleanup();
  });
  header.querySelector(".variant-close")?.addEventListener("click", cleanup);
  document.addEventListener("keydown", escHandler);

  return { body, close: cleanup };
}

function formatPredictionCell(tool, rawValue, store, meta = {}) {
  if (!tool) return displayValue(rawValue);
  if (String(tool).toLowerCase() !== "expectosc" || !rawValue || rawValue === "—") {
    return displayValue(rawValue);
  }
  const parsed = typeof rawValue === "string" ? safeJson(rawValue) : rawValue;
  if (!parsed) return displayValue(rawValue);
  const scores = Array.isArray(parsed.scores)
    ? parsed.scores
        .map((rec) => ({
          tissue: rec?.tissue || "",
          score: Number(rec?.score),
        }))
        .filter((rec) => rec.tissue && Number.isFinite(rec.score))
    : [];
  const topUp = [...scores].filter((s) => s.score > 0).sort((a, b) => b.score - a.score)[0] || null;
  const topDown = [...scores].filter((s) => s.score < 0).sort((a, b) => a.score - b.score)[0] || null;
  const id = `exp_${rand7()}`;
  store.set(id, {
    topUp,
    topDown,
    topList: Array.isArray(parsed.top_tissues) ? parsed.top_tissues : [],
    bottomList: Array.isArray(parsed.bottom_tissues) ? parsed.bottom_tissues : [],
    scores,
    variant_id: meta.variantId || "",
    variant_label: meta.variantLabel || meta.variantId || "",
    location_label: meta.locationLabel || "",
    gene_label: meta.geneLabel || "",
  });
  return {
    __html: renderExpectoCellHtml(id, topUp, topDown),
  };
}

function renderExpectoCellHtml(id, topUp, topDown) {
  const upScore = topUp ? topUp.score : null;
  const downScore = topDown ? topDown.score : null;
  const upColor = expectoChipColor(upScore);
  const downColor = expectoChipColor(downScore);
  const upText = topUp ? `↑ ${formatScore(upScore)}` : "↑ —";
  const downText = topDown ? `↓ ${formatScore(downScore)}` : "↓ —";
  return `
    <div class="expecto-cell" data-exp-id="${escapeInline(id)}" style="display:flex;gap:8px;">
      <span class="exp-chip exp-up" data-exp-toggle style="flex:1;display:flex;justify-content:center;padding:4px 10px;border-radius:999px;font-size:12px;font-weight:600;cursor:pointer;background:${upColor.bg};color:${upColor.fg};border:${upColor.border};">
        ${escapeInline(upText)}
      </span>
      <span class="exp-chip exp-down" data-exp-toggle style="flex:1;display:flex;justify-content:center;padding:4px 10px;border-radius:999px;font-size:12px;font-weight:600;cursor:pointer;background:${downColor.bg};color:${downColor.fg};border:${downColor.border};">
        ${escapeInline(downText)}
      </span>
    </div>
  `;
}

function wireExpectoPanels(root, store) {
  if (!store || store.size === 0) return;
  root.querySelectorAll(".expecto-cell").forEach((cell) => {
    const id = cell.getAttribute("data-exp-id");
    const data = store.get(id);
    if (!id || !data) return;
    cell.querySelectorAll("[data-exp-toggle]").forEach((btn) => {
      btn.addEventListener("click", () => openExpectoModal(data));
    });
  });
}

function renderExpectoPanel(data) {
  const lists = `
    <div style="display:flex;gap:16px;flex-wrap:wrap;">
      ${renderExpectoList("Top up-regulated cell types", data.topList)}
      ${renderExpectoList("Top down-regulated cell types", data.bottomList)}
    </div>
  `;
  return `
    <div class="expecto-panel-inner" style="padding:8px 10px;border-radius:12px;background:rgba(255,255,255,0.03);border:1px solid var(--line);">
      ${lists}
      <div style="margin-top:10px;">
        <div class="muted" style="font-size:12px;letter-spacing:0.5px;text-transform:uppercase;">Effect across cell types</div>
        <div class="exp-bars-controls" style="margin:8px 0;display:flex;gap:8px;flex-wrap:wrap;">
          <button type="button" class="exp-bars-btn active" data-exp-filter="gte1" style="padding:4px 12px;border-radius:999px;border:1px solid var(--line);background:rgba(148,163,184,0.18);color:var(--fg);cursor:pointer;font-size:12px;font-weight:600;">|score| ≥ 1</button>
          <button type="button" class="exp-bars-btn" data-exp-filter="all" style="padding:4px 12px;border-radius:999px;border:1px solid var(--line);background:transparent;color:var(--fg);cursor:pointer;font-size:12px;font-weight:600;">Show all cell types</button>
        </div>
        <div class="exp-bars" data-exp-scores='${jsonAttr(data.scores || [])}' data-current-filter="gte1" data-variant-label="${escapeInline(data.variant_label || data.variant_id || "")}" data-location-label="${escapeInline(data.location_label || "")}" data-gene-label="${escapeInline(data.gene_label || "")}" style="margin-top:6px;"></div>
        <div class="exp-bars-actions" style="margin-top:10px;display:flex;gap:8px;flex-wrap:wrap;">
          <button type="button" class="exp-save-btn" data-save-mode="current" style="padding:4px 12px;border-radius:999px;border:1px solid var(--line);background:rgba(148,163,184,0.12);color:var(--fg);cursor:pointer;font-size:12px;font-weight:600;">Save current view</button>
          <button type="button" class="exp-save-btn" data-save-mode="all" style="padding:4px 12px;border-radius:999px;border:1px solid var(--line);background:rgba(148,163,184,0.12);color:var(--fg);cursor:pointer;font-size:12px;font-weight:600;">Save all cell types</button>
        </div>
      </div>
    </div>
  `;
}

function renderExpectoList(title, items) {
  const safeItems = (items || []).filter(Boolean).slice(0, 4);
  if (!safeItems.length) return "";
  const rows = safeItems
    .map((rec) => {
      const tissue = formatCellTissueLabel(rec?.tissue || rec?.cell_type || "");
      const score = Number(rec?.score);
      return `
        <div style="display:flex;align-items:center;gap:6px;">
          <span style="font-weight:600;">${scoreLabel(score)}</span>
          <span>${escapeInline(String(tissue || "Cell type not provided"))}</span>
        </div>
      `;
    })
    .join("");
  return `
    <div style="flex:1;min-width:220px;">
      <div class="muted" style="font-size:12px;letter-spacing:0.5px;text-transform:uppercase;">${escapeInline(title)}</div>
      <div style="margin-top:6px;display:flex;flex-direction:column;gap:4px;">${rows}</div>
    </div>
  `;
}

function renderExpectoBars(scores, mode = "gte1") {
  const dataset = prepareExpectoDataset(scores, mode);
  if (!dataset.rows.length) {
    return mode === "gte1"
      ? `<div class="muted">No cell types with |score| ≥ 1. Switch to "Show all".</div>`
      : `<div class="muted">No ExpectoSC score distribution.</div>`;
  }
  const rows = dataset.rows
    .map((rec) => {
      const tissue = rec.tissue || "";
      const displayLabel = formatCellTissueLabel(rec.tissue || "");
      const score = rec.score;
      const pct = Math.min(48, (Math.abs(score) / dataset.maxAbs) * 48);
      const left = score >= 0 ? 52 : 52 - pct;
      const color = Number.isFinite(score) && Math.abs(score) >= 1 ? tissueColor(tissue) : "var(--muted)";
      const bg = Number.isFinite(score) && Math.abs(score) >= 1 ? color : "var(--muted)";
      return `
        <div class="exp-bar-row" style="display:flex;align-items:center;gap:10px;font-size:12px;padding:3px 0;">
          <div style="flex:0 0 220px;">${escapeInline(displayLabel)}</div>
          <div style="flex:1;position:relative;height:14px;background:rgba(148,163,184,0.15);border-radius:8px;overflow:hidden;">
            <div style="position:absolute;top:0;bottom:0;width:2px;background:rgba(255,255,255,0.5);left:52%;transform:translateX(-1px);opacity:0.6;"></div>
            <div style="position:absolute;top:2px;bottom:2px;border-radius:8px;background:${Number.isFinite(score) && Math.abs(score) >= 1 ? color : "rgba(148,163,184,0.6)"};left:${left}%;width:${pct}%;"></div>
          </div>
          <div style="width:60px;text-align:right;font-variant-numeric:tabular-nums;color:${Number.isFinite(score) && Math.abs(score) >= 1 ? color : "var(--muted)"};">${formatScore(score)}</div>
        </div>
      `;
    })
    .join("");
  return `
    <div style="position:relative;padding:8px 0;">
      <div style="position:absolute;top:0;bottom:0;left:52%;width:1px;background:rgba(255,255,255,0.2);"></div>
      <div style="display:flex;flex-direction:column;gap:2px;max-height:260px;overflow:auto;padding-right:4px;">
        ${rows}
      </div>
      <div style="display:flex;justify-content:space-between;font-size:11px;color:var(--muted);margin-top:4px;">
        <span>Negative</span>
        <span>Positive</span>
      </div>
      ${renderTissueLegend(dataset.rows)}
    </div>
  `;
}

function scoreLabel(score) {
  if (!Number.isFinite(score)) return "—";
  return score >= 0 ? `↑ ${formatScore(score)}` : `↓ ${formatScore(score)}`;
}

function renderMiniTable(header, rows, columnStyles = []) {
  const styleAttr = (idx) => {
    const style = columnStyles[idx];
    return style ? ` style="${defaultEscapeHtml(String(style))}"` : "";
  };
  const thead = `<thead><tr>${header.map((h, idx) => `<th class="wrap"${styleAttr(idx)}>${escapeInline(h)}</th>`).join("")}</tr></thead>`;
  const tbody = `<tbody>${rows.map((r) => `<tr>${r.map((v, idx) => `<td class="wrap"${styleAttr(idx)}>${renderCellValue(v)}</td>`).join("")}</tr>`).join("")}</tbody>`;
  return `<div style="overflow:auto; max-height:420px;"><table class="vs-table" style="table-layout:auto;width:100%;">${thead}${tbody}</table></div>`;
}

function renderCellValue(v) {
  if (v && typeof v === "object" && v.__html) return v.__html;
  return escapeInline(v);
}

function safeJson(str) {
  if (typeof str !== "string") return null;
  const unescaped = str
    .replace(/&quot;/g, '"')
    .replace(/&#39;/g, "'")
    .replace(/&lt;/g, "<")
    .replace(/&gt;/g, ">")
    .replace(/&amp;/g, "&");
  try {
    return JSON.parse(unescaped);
  } catch {
    return null;
  }
}

function formatScore(score) {
  if (!Number.isFinite(score)) return "—";
  return Number(score).toFixed(2);
}

function expectoChipColor(score) {
  if (!Number.isFinite(score)) {
    return { bg: "rgba(148,163,184,0.25)", fg: "var(--fg)", border: "1px solid rgba(148,163,184,0.4)" };
  }
  if (Math.abs(score) < 1) {
    return { bg: "rgba(148,163,184,0.2)", fg: "var(--fg)", border: "1px solid rgba(148,163,184,0.5)" };
  }
  if (score >= 0) {
    return { bg: "rgba(34,197,94,0.2)", fg: "var(--fg)", border: "1px solid rgba(34,197,94,0.45)" };
  }
  return { bg: "rgba(248,113,113,0.2)", fg: "var(--fg)", border: "1px solid rgba(248,113,113,0.45)" };
}

function openExpectoModal(data) {
  closeExpectoModal();
  const overlay = document.createElement("div");
  overlay.className = "expecto-overlay";
  overlay.style.position = "fixed";
  overlay.style.inset = "0";
  overlay.style.background = "rgba(2,4,12,0.65)";
  overlay.style.backdropFilter = "blur(4px)";
  overlay.style.zIndex = "9999";
  overlay.style.display = "flex";
  overlay.style.alignItems = "center";
  overlay.style.justifyContent = "center";

  const modal = document.createElement("div");
  modal.className = "expecto-modal";
  modal.style.maxWidth = "760px";
  modal.style.width = "90%";
  modal.style.maxHeight = "85vh";
  modal.style.overflow = "auto";
  modal.style.background = "var(--card)";
  modal.style.border = "1px solid var(--line)";
  modal.style.borderRadius = "18px";
  modal.style.boxShadow = "0 30px 80px rgba(0,0,0,0.45)";
  modal.style.padding = "18px 20px";
  modal.innerHTML = `
    <div style="display:flex;align-items:center;justify-content:space-between;gap:16px;">
      <div>
        <div style="font-weight:700;">ExpectoSC cell type effects</div>
        <div class="muted" style="font-size:13px;">
          ${escapeInline(data.variant_label || data.variant_id || "")}
          ${data.location_label ? ` • ${escapeInline(data.location_label)}` : ""}
          ${data.gene_label ? ` • ${escapeInline(data.gene_label)}` : ""}
        </div>
      </div>
      <button type="button" class="expecto-close" style="border:none;background:transparent;font-size:20px;cursor:pointer;color:var(--fg);">×</button>
    </div>
    <div style="margin-top:12px;">${renderExpectoPanel(data)}</div>
  `;

  overlay.appendChild(modal);
  document.body.appendChild(overlay);

  const closeBtn = modal.querySelector(".expecto-close");
  const cleanup = () => {
    overlay.remove();
    document.removeEventListener("keydown", escHandler);
  };
  const escHandler = (e) => {
    if (e.key === "Escape") cleanup();
  };
  overlay.addEventListener("click", (e) => {
    if (e.target === overlay) cleanup();
  });
  closeBtn?.addEventListener("click", cleanup);
  document.addEventListener("keydown", escHandler);
  overlay.setAttribute("data-exp-overlay", "1");
  initExpectoModal(modal);
}

function closeExpectoModal() {
  document.querySelectorAll('[data-exp-overlay="1"]').forEach((el) => el.remove());
}

function initExpectoModal(modal) {
  modal.querySelectorAll(".exp-bars").forEach((container) => {
    const scores = safeJson(container.getAttribute("data-exp-scores")) || [];
    const meta = {
      variantLabel: container.getAttribute("data-variant-label") || "",
      locationLabel: container.getAttribute("data-location-label") || "",
      geneLabel: container.getAttribute("data-gene-label") || "",
    };
    const render = (mode) => {
      container.innerHTML = renderExpectoBars(scores, mode);
      container.setAttribute("data-current-filter", mode);
    };
    render("gte1");
    const controls = container.parentElement?.querySelectorAll(".exp-bars-btn") || [];
    controls.forEach((btn) => {
      btn.addEventListener("click", () => {
        const mode = btn.getAttribute("data-exp-filter") || "gte1";
        controls.forEach((b) => setExpectoFilterButtonState(b, b === btn));
        render(mode);
      });
      setExpectoFilterButtonState(btn, btn.classList.contains("active"));
    });

    const saveButtons = container.parentElement?.querySelectorAll(".exp-save-btn") || [];
    saveButtons.forEach((btn) => {
      btn.addEventListener("click", () => {
        const mode = btn.getAttribute("data-save-mode") === "all" ? "all" : (container.getAttribute("data-current-filter") || "gte1");
        exportExpectoChart(scores, mode, meta);
      });
    });
  });
}

function renderTissueLegend(rows) {
  const seen = new Set();
  const items = [];
  for (const rec of rows) {
    const tissue = extractTissuePart(rec.tissue || "");
    if (!tissue || seen.has(tissue)) continue;
    seen.add(tissue);
    const color = Math.abs(rec.score) >= 1 ? tissueColor(tissue) : "var(--muted)";
    items.push({ label: tissue, color });
    if (items.length >= 6) break;
  }
  if (!items.length) return "";
  const legendRows = items
    .map((item) => `
      <div style="display:flex;align-items:center;gap:6px;font-size:12px;color:var(--muted);">
        <span style="width:10px;height:10px;border-radius:999px;background:${item.color};display:inline-block;"></span>
        <span>${escapeInline(item.label)}</span>
      </div>
    `)
    .join("");
  return `
    <div style="margin-top:10px;">
      <div class="muted" style="font-size:11px;text-transform:uppercase;letter-spacing:0.5px;">Color legend (tissue)</div>
      <div style="margin-top:4px;display:grid;grid-template-columns:repeat(auto-fit,minmax(180px,1fr));gap:4px;">
        ${legendRows}
      </div>
    </div>
  `;
}

function setExpectoFilterButtonState(btn, active) {
  if (!btn) return;
  btn.classList.toggle("active", !!active);
  if (active) {
    btn.style.background = "rgba(148,163,184,0.22)";
    btn.style.color = "var(--fg)";
    btn.style.borderColor = "rgba(148,163,184,0.6)";
  } else {
    btn.style.background = "rgba(148,163,184,0.08)";
    btn.style.color = "var(--muted)";
    btn.style.borderColor = "rgba(148,163,184,0.3)";
  }
}

function prepareExpectoDataset(scores, mode = "gte1") {
  const arr = Array.isArray(scores)
    ? scores.map((rec) => ({
        tissue: rec?.tissue || rec?.cell_type || "",
        score: Number(rec?.score),
      })).filter((rec) => rec.tissue && Number.isFinite(rec.score))
    : [];
  if (!arr.length) return { rows: [], maxAbs: 1 };
  const filtered = mode === "gte1" ? arr.filter((r) => Math.abs(r.score) >= 1) : arr;
  if (!filtered.length) return { rows: [], maxAbs: 1 };
  const sorted = [...filtered].sort((a, b) => Math.abs(b.score) - Math.abs(a.score));
  const limited = mode === "all" ? sorted : sorted.slice(0, 30);
  const maxAbs = Math.max(...sorted.map((r) => Math.abs(r.score))) || 1;
  return { rows: limited, maxAbs };
}

function exportExpectoChart(scores, mode, meta) {
  const dataset = prepareExpectoDataset(scores, mode);
  if (!dataset.rows.length) {
    alert("No data to export for this view.");
    return;
  }
  const legendItems = collectLegendItems(dataset.rows);
  const svg = buildExpectoSvg(dataset.rows, dataset.maxAbs, meta, legendItems, mode);
  const name = `${(meta.variantLabel || "variant").replace(/\s+/g, "_")}_${mode}_expectosc.svg`;
  downloadText(svg, name, "image/svg+xml");
}

function collectLegendItems(rows) {
  const seen = new Set();
  const items = [];
  for (const rec of rows) {
    const tissue = extractTissuePart(rec.tissue || "");
    if (!tissue || seen.has(tissue)) continue;
    seen.add(tissue);
    items.push({ label: tissue, color: Math.abs(rec.score) >= 1 ? tissueColor(tissue) : "var(--muted)" });
    if (items.length >= 10) break;
  }
  return items;
}

function buildExpectoSvg(rows, maxAbs, meta, legendItems, mode) {
  const width = 900;
  const height = Math.max(200, rows.length * 26 + 200);
  const chartWidth = 520;
  const chartHeight = rows.length * 24;
  const marginLeft = 230;
  const zeroX = marginLeft + chartWidth / 2;

  const title = `${meta.variantLabel || ""}${meta.locationLabel ? ` • ${meta.locationLabel}` : ""}${meta.geneLabel ? ` • ${meta.geneLabel}` : ""}`.trim();
  const subtitle = mode === "all" ? "All cell types" : "Cell types with |score| ≥ 1";

  const bars = rows
    .map((rec, idx) => {
      const y = 100 + idx * 24;
      const tissue = rec.tissue || "";
      const label = formatCellTissueLabel(tissue);
      const score = rec.score;
      const barWidth = (Math.abs(score) / maxAbs) * (chartWidth / 2);
      const x = score >= 0 ? zeroX : zeroX - barWidth;
      const color = Math.abs(score) >= 1 ? tissueColor(tissue) : "#94a3b8";
      return `
        <text x="${marginLeft - 10}" y="${y + 12}" font-size="12" text-anchor="end" fill="#cbd5f5">${escapeSvg(label)}</text>
        <rect x="${x}" y="${y}" width="${barWidth}" height="14" rx="7" fill="${color}"></rect>
        <text x="${score >= 0 ? x + barWidth + 6 : x - 6}" y="${y + 12}" font-size="12" text-anchor="${score >= 0 ? "start" : "end"}" fill="#e2e8f0">${formatScore(score)}</text>
      `;
    })
    .join("");

  const legend = legendItems.length
    ? legendItems
        .map((item, i) => `
          <g transform="translate(${marginLeft}, ${chartHeight + 140 + i * 18})">
            <rect width="10" height="10" fill="${item.color}"></rect>
            <text x="16" y="10" font-size="12" fill="#cbd5f5">${escapeSvg(item.label)}</text>
          </g>
        `)
        .join("")
    : "";

  return `
    <svg xmlns="http://www.w3.org/2000/svg" width="${width}" height="${height}" viewBox="0 0 ${width} ${height}" style="background:#0b1020;font-family:'Inter', system-ui, sans-serif;">
      <text x="${marginLeft}" y="40" font-size="20" font-weight="700" fill="#f8fafc">${escapeSvg(title || "ExpectoSC cell type effects")}</text>
      <text x="${marginLeft}" y="60" font-size="14" fill="#cbd5f5">${escapeSvg(subtitle)}</text>
      <line x1="${zeroX}" y1="90" x2="${zeroX}" y2="${chartHeight + 110}" stroke="rgba(255,255,255,0.3)" stroke-dasharray="4 3"></line>
      ${bars}
      <text x="${marginLeft}" y="${chartHeight + 120}" font-size="12" fill="#cbd5f5">Negative</text>
      <text x="${marginLeft + chartWidth}" y="${chartHeight + 120}" font-size="12" text-anchor="end" fill="#cbd5f5">Positive</text>
      ${legend}
    </svg>
  `;
}

function downloadText(content, filename, mime) {
  const blob = new Blob([content], { type: mime || "text/plain" });
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url;
  a.download = filename;
  document.body.appendChild(a);
  a.click();
  document.body.removeChild(a);
  setTimeout(() => URL.revokeObjectURL(url), 500);
}

function tissueColor(label) {
  const tissue = extractTissuePart(label);
  const hue = hashString(tissue) % 360;
  return `hsl(${hue}, 65%, 55%)`;
}

function extractTissuePart(label) {
  const parts = String(label || "").split("/");
  return (parts[1] || parts[0] || "").trim() || "tissue";
}

function splitCellTissue(label) {
  const parts = String(label || "").split("/");
  return {
    cell: (parts[0] || "").trim(),
    tissue: (parts[1] || "").trim(),
  };
}

function formatCellTissueLabel(label) {
  const { cell, tissue } = splitCellTissue(label);
  if (cell && tissue) return `${cell} (${tissue})`;
  return cell || tissue || "";
}

function hashString(str) {
  let hash = 0;
  for (let i = 0; i < str.length; i++) {
    hash = (hash << 5) - hash + str.charCodeAt(i);
    hash |= 0;
  }
  return Math.abs(hash);
}

function rand7() {
  return Math.random().toString(36).slice(2, 9);
}

function jsonAttr(obj) {
  try {
    return defaultEscapeHtml(JSON.stringify(obj));
  } catch {
    return "[]";
  }
}

function defaultEscapeHtml(s) {
  return String(s ?? "").replace(/[&<>"']/g, (m) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[m]));
}
function escapeInline(s) {
  return defaultEscapeHtml(String(s ?? ""));
}

function escapeSvg(s) {
  return String(s ?? "").replace(/[&<>"']/g, (m) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[m]));
}

function defaultParseCSV(text, delimiter = ",") {
  if (!text) return [];
  const lines = String(text).replace(/\r\n?/g, "\n").split("\n").filter((l) => l.trim());
  if (!lines.length) return [];
  const parseLine = (line) => {
    const out = [];
    let cur = "";
    let inQ = false;
    for (let i = 0; i < line.length; i++) {
      const ch = line[i];
      if (inQ) {
        if (ch === '"' && line[i + 1] === '"') {
          cur += '"';
          i++;
        } else if (ch === '"') {
          inQ = false;
        } else {
          cur += ch;
        }
      } else {
        if (ch === '"') inQ = true;
        else if (ch === delimiter) {
          out.push(cur);
          cur = "";
        } else {
          cur += ch;
        }
      }
    }
    out.push(cur);
    return out;
  };
  const header = parseLine(lines[0]).map((h) => h.trim());
  return lines.slice(1).map((raw) => {
    const parts = parseLine(raw);
    const obj = {};
    for (let i = 0; i < header.length; i++) obj[header[i]] = (parts[i] ?? "").trim();
    return obj;
  });
}

function uniquePreserve(items) {
  const seen = new Set();
  const out = [];
  for (const item of items) {
    if (!seen.has(item)) {
      seen.add(item);
      out.push(item);
    }
  }
  return out;
}

function normVariant(s) {
  const v = (s || "").trim();
  return v || "";
}
function normGene(s) {
  const v = (s || "").trim();
  return v ? v.toUpperCase() : "";
}

function formatChrom(chrom) {
  const c = String(chrom || "").trim();
  return c ? (c.startsWith("chr") ? c : `chr${c}`) : "—";
}
function formatPos(pos) {
  if (pos == null || pos === "") return "—";
  const num = Number(pos);
  if (Number.isFinite(num)) return num.toLocaleString("en-US");
  return String(pos);
}

function formatLocationLabel(chrom, pos) {
  const chromStr = formatChrom(chrom);
  const posStr = pos == null || pos === "" ? "" : String(pos).replace(/,/g, "");
  if (!posStr || chromStr === "—") return "—";
  return `${chromStr}:${posStr}`;
}

function displayValue(val) {
  if (val == null) return "—";
  const str = String(val).trim();
  if (!str) return "—";
  return str.replace(/_/g, " ");
}

function toolLabel(name) {
  return name ? name.replace(/_/g, " ") : "";
}

function parseNumber(val) {
  if (val == null || val === "") return null;
  const num = Number(val);
  return Number.isFinite(num) ? num : null;
}


function dedupeAlleleRows(rows) {
  const seen = new Set();
  const result = [];
  for (const rec of rows) {
    const key = `${rec.study || ""}|${rec.ref || ""}|${rec.alt || ""}`;
    if (seen.has(key)) continue;
    seen.add(key);
    result.push(rec);
  }
  return result;
}

function renderAlleleChange(ref, alt) {
  if (!ref && !alt) return "—";
  if (!ref) return `→ ${alt}`;
  if (!alt) return `${ref} →`;
  return `${ref}→${alt}`;
}

function sortTraits(traits) {
  return [...traits].sort((a, b) => {
    const riskA = parseNumber(a.risk_score);
    const riskB = parseNumber(b.risk_score);
    if (Number.isFinite(riskA) || Number.isFinite(riskB)) {
      return (riskB ?? -Infinity) - (riskA ?? -Infinity);
    }
    const pA = parseNumber(a.p_value);
    const pB = parseNumber(b.p_value);
    if (Number.isFinite(pA) || Number.isFinite(pB)) {
      return (pA ?? Infinity) - (pB ?? Infinity);
    }
    return 0;
  });
}

function formatPercent(val) {
  if (!Number.isFinite(val)) return "—";
  return `${(val * 100).toFixed(2).replace(/\.00$/, "")}%`;
}

function formatInteger(val) {
  if (!Number.isFinite(val)) return "—";
  return Number(val).toLocaleString("en-US");
}

function formatRisk(val) {
  if (!Number.isFinite(val)) return "—";
  if (Math.abs(val) >= 100) return val.toFixed(1);
  return val.toFixed(3).replace(/0+$/, "").replace(/\.$/, "");
}

function formatPValue(val) {
  if (!Number.isFinite(val)) return "—";
  if (val === 0) return "0";
  if (val < 1e-4) return val.toExponential(2);
  return val.toPrecision(3);
}
