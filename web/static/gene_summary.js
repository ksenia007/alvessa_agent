// Ksenia Sokolova
// 2025 - 10 - 22
// static/gene_summary.js
// "Gene Summary" card.

// Expects directory structure (case-sensitive symbols):
//   genes/genes.index.csv
//   genes/<SYMBOL>/gene.json
//   genes/<SYMBOL>/transcripts.tsv
//   genes/<SYMBOL>/interactions_human.tsv (optional)
//   genes/<SYMBOL>/interactions_nonhuman.tsv (optional)
//   genes/<SYMBOL>/summary.txt (optional)


const GENE_RENDER_MODE = "chips"; // "chips" | "comma"
const SECT_SUMMARY_STYLE = "class='muted sect-hdr sub'";
const TOP_SUMMARY_STYLE = "class='sect-hdr top-hdr'";  


function renderGeneList(genes) {
  if (!Array.isArray(genes) || genes.length === 0) return "";
  if (GENE_RENDER_MODE === "comma") {
    return defaultEscapeHtml(genes.join(", "));
  }
  // chips
  return genes.map(g => `<span class="chip wrap">${defaultEscapeHtml(String(g))}</span>`).join("");
}


function msigdbCategoryNote(cat) {
    const m = {
      H:  "Hallmark gene sets (broad biological states/processes)",
      C1: "Positional gene sets (cytogenetic bands)",
      C2: "Curated gene sets (pathways/literature)",
      C3: "Regulatory target sets (TF & miRNA targets)",
      C4: "Computational gene sets (co-expression/other)",
      C5: "Gene Ontology (GO) sets",
      C6: "Oncogenic signatures",
      C7: "Immunologic signatures",
      C8: "Cell type signatures (single-cell)",
    };
    return m[cat] || "";
  }
  
  function renderMsigdbBlock(msigdb_annotations) {
    if (!msigdb_annotations || typeof msigdb_annotations !== "object") return "";
    const order = ["H","C1","C2","C3","C4","C5","C6","C7","C8"];
    const inner = [];
  
    for (const cat of order) {
      const arr = msigdb_annotations[cat];
      if (!Array.isArray(arr) || arr.length === 0) continue;
  
      const list = arr.map(x => defaultEscapeHtml(String(x))).join(", ");
      const note = msigdbCategoryNote(cat);
  
      inner.push(`
        <details style="margin-top:6px;">
          <summary ${SECT_SUMMARY_STYLE}>
            ${`MSigDB ${cat}${note ? ` — ${defaultEscapeHtml(note)}` : ""}`} (${arr.length})
          </summary>
          <div style="margin-top:6px; white-space:normal; line-height:1.4;">
            ${list}
          </div>
        </details>
      `);
    }
  
    if (!inner.length) return "";
    return `
      <details style="margin-top:8px;">
        <summary ${TOP_SUMMARY_STYLE}>MSigDB</summary>
        <div style="margin-top:6px;">
          ${inner.join("\n")}
        </div>
      </details>
    `;
  }
  
  
  function isNonEmptyObject(o){
    return o && typeof o === "object" && Object.keys(o).length > 0;
  }
  
  function normalizeDiseaseAnnotations(x){
    // Accept array OR dict {'-': 'A', '1': 'B', ...} -> return array of values
    if (!x) return [];
    if (Array.isArray(x)) return x.filter(Boolean).map(v => String(v));
    if (typeof x === "object") {
      return Object.values(x).filter(Boolean).map(v => String(v));
    }
    return [];
  }

  function renderTranscriptsSection(obj) {
    if (!obj) return "";
    // Accept map {ENST: {...}} or array of entries
    let rows = [];
  
    if (Array.isArray(obj)) {
      // Array could be ["ENST...", {id:"ENST...", n_exons:...}, ...]
      for (const it of obj) {
        if (!it) continue;
        if (typeof it === "string") {
          rows.push({ id: it, n_exons: "", other: "" });
        } else if (typeof it === "object") {
          const id = it.id || it.transcript_id || it.name || "";
          const { n_exons, ...rest } = it;
          rows.push({
            id: String(id || "").trim(),
            n_exons: (n_exons ?? "") === "" ? "" : String(n_exons),
            other: summarizeRest(rest)
          });
        }
      }
    } else if (typeof obj === "object") {
      // Map of id -> { n_exons, ... }
      for (const [id, v] of Object.entries(obj)) {
        if (!v || typeof v !== "object") {
          rows.push({ id, n_exons: "", other: "" });
          continue;
        }
        const { n_exons, ...rest } = v;
        rows.push({
          id,
          n_exons: (n_exons ?? "") === "" ? "" : String(n_exons),
          other: summarizeRest(rest)
        });
      }
    }
  
    if (!rows.length) return "";
  
    // Sort by ID asc
    rows.sort((a, b) => a.id.localeCompare(b.id));
  
    const header = ["Transcript ID", "Exons", "Other"];
    const body = rows.map(r => [
      escapeInline(r.id),
      escapeInline(r.n_exons ?? ""),
      escapeInline(r.other ?? "")
    ]);
  
    return `
      <details style="margin-top:8px;">
        <summary ${TOP_SUMMARY_STYLE}>Transcripts (${rows.length})</summary>
        <div style="margin-top:6px;">
          ${renderMiniTable(header, body)}
        </div>
      </details>
    `;
  
    function summarizeRest(rest) {
      // Show compact K:V; arrays joined, nested objects shown as JSON once
      const parts = [];
      for (const [k, v] of Object.entries(rest)) {
        if (v == null) continue;
        if (Array.isArray(v)) parts.push(`${k}: ${v.join(", ")}`);
        else if (typeof v === "object") parts.push(`${k}: ${JSON.stringify(v)}`);
        else parts.push(`${k}: ${String(v)}`);
      }
      return parts.join(" • ");
    }
  }
  
  function renderIsoformsSection(obj) {
    if (!obj || typeof obj !== "object") return "";
  
    // 1) Handle general_localization (either key or inline object)
    const gl =
      obj.general_localization ||
      (obj.name === "general_localization" ? obj : null);
  
    const glBlock = gl ? renderGeneralLocalizationSection("Isoforms", { general_localization: gl }) : "";
  
    // 2) Collect isoform entries: keys like "P38398-3" etc (exclude general_localization)
    const entries = [];
    for (const [key, v] of Object.entries(obj)) {
      if (key === "general_localization") continue;
      if (!v || typeof v !== "object") continue;
      entries.push([key, v]);
    }
    if (!glBlock && entries.length === 0) return "";
  
    entries.sort((a, b) => a[0].localeCompare(b[0]));
  
    const cards = entries.map(([id, info]) => {
      const name    = info.name ? ` — ${escapeInline(info.name)}` : "";
      const status  = info.status ? `<span class="chip">Status: ${escapeInline(info.status)}</span>` : "";
      const aliases = Array.isArray(info.aliases) && info.aliases.length
        ? `<div class="row" style="gap:6px;flex-wrap:wrap;margin-top:4px;">
             ${info.aliases.map(a => `<span class="chip">${escapeInline(String(a))}</span>`).join("")}
           </div>`
        : "";
      const locs = Array.isArray(info.locations) && info.locations.length
        ? `<div class="row" style="gap:6px;flex-wrap:wrap;margin-top:6px;">
             ${info.locations.map(a => `<span class="chip">${escapeInline(String(a))}</span>`).join("")}
           </div>`
        : "";
      const notes = Array.isArray(info.notes) && info.notes.length
        ? `<div class="row" style="gap:6px;flex-wrap:wrap;margin-top:6px;">
             ${info.notes.map(n => `<span class="chip chip-note">${escapeInline(String(n))}</span>`).join("")}
           </div>`
        : "";
      return `
        <div class="card" style="margin:8px 0; padding:12px;">
          <div class="muted" style="font-weight:700;">${escapeInline(id)}${name}</div>
          ${status ? `<div class="row" style="gap:6px;flex-wrap:wrap;margin-top:4px;">${status}</div>` : ""}
          ${aliases}
          ${locs}
          ${notes}
        </div>
      `;
    }).join("");
  
    const total = entries.length;
  
    return `
      <details style="margin-top:8px;">
        <summary ${TOP_SUMMARY_STYLE}>Isoforms${total ? ` (${total})` : ""}</summary>
        <div style="margin-top:6px;">
          ${glBlock}
          ${cards}
        </div>
      </details>
    `;
  }
  
  function renderGeneralLocalizationSection(title, obj){
    const gl = obj.general_localization;
    if (!gl) return "";
    const locs  = Array.isArray(gl.locations) ? gl.locations : [];
    const notes = Array.isArray(gl.notes) ? gl.notes : [];
  
    const locRow  = locs.length
      ? `<div class="row" style="gap:6px;flex-wrap:wrap;">${locs.map(x=>`<span class="chip">${escapeInline(String(x))}</span>`).join("")}</div>`
      : "";
  
    const noteRow = notes.length
      ? `<div class="row" style="gap:6px;flex-wrap:wrap;margin-top:6px;">${notes.map(n=>`<span class="chip chip-note">${escapeInline(String(n))}</span>`).join("")}</div>`
      : "";
  
    const meta = [];
    if (gl.status) meta.push(`Status: <strong>${escapeInline(String(gl.status))}</strong>`);
    if (Array.isArray(gl.aliases) && gl.aliases.length) meta.push(`Aliases: ${escapeInline(gl.aliases.join(", "))}`);
  
    return `
      <div class="card" style="margin:8px 0; padding:12px;">
        <div class="muted" style="font-weight:700; margin-bottom:6px;">${escapeInline(title)} - General localization</div>
        ${meta.length ? `<div class="muted" style="margin-bottom:6px;">${meta.join(" • ")}</div>` : ""}
        ${locRow || `<div class="muted">No locations reported.</div>`}
        ${noteRow}
      </div>
    `;
  }
  
  
  function renderOpenTargetsBlock(ot){
    if (!ot || typeof ot !== "object") return "";
  
    const diseases = normalizeDiseaseAnnotations(ot.disease_annotations);
    const diseaseInner = diseases.length ? `
      <details style="margin-top:6px;">
        <summary ${SECT_SUMMARY_STYLE}>Disease annotations (${diseases.length})</summary>
        <div style="margin-top:6px; white-space:normal; line-height:1.4;">
          ${defaultEscapeHtml(diseases.join(", "))}
        </div>
      </details>` : "";
  
    const otherBlocks = [];
    for (const [k, v] of Object.entries(ot)) {
      if (k === "is_essential" || k === "disease_annotations") continue;
      const hasVal = Array.isArray(v) ? v.length > 0 : isNonEmptyObject(v);
      if (!hasVal) continue;
  
      if (Array.isArray(v)) {
        const list = v.map(x => `<span class="pill wrap">${defaultEscapeHtml(String(x))}</span>`).join(" ");
        otherBlocks.push(`
          <details style="margin-top:6px;">
            <summary ${SECT_SUMMARY_STYLE}>${defaultEscapeHtml(k)}</summary>
            <div style="margin-top:6px;">${list}</div>
          </details>
        `);
      } else if (typeof v === "object") {
        otherBlocks.push(`
          <details style="margin-top:6px;">
            <summary ${SECT_SUMMARY_STYLE}>${defaultEscapeHtml(k)}</summary>
            <div style="margin-top:6px;">${objectToTable(v, defaultEscapeHtml)}</div>
          </details>
        `);
      }
    }
  
    const essentialBadge = (typeof ot.is_essential === "boolean")
    ? `<span class="pill wrap" title="from Open Targets">is_essential:&nbsp;<strong>${ot.is_essential ? "true" : "false"}</strong></span>`    
      : "";
  
    const inner = [diseaseInner, ...otherBlocks].filter(Boolean).join("\n");
    if (!inner) return "";
  
    return `
      <details style="margin-top:8px;">
        <summary ${TOP_SUMMARY_STYLE}>Open Targets</summary>
        <div style="margin-top:6px;">
        ${essentialBadge ? `<div class="row" style="gap:6px;flex-wrap:wrap;margin-bottom:6px;">${essentialBadge}</div>` : ""}
           ${inner}
        </div>
      </details>
    `;
  }

  function renderGwasTraitsBlock(gwasObj, traitsArr){
    const hasGwas = gwasObj && typeof gwasObj === "object" && Object.keys(gwasObj).length;
    const hasTraits = Array.isArray(traitsArr) && traitsArr.length > 0;
  
    if (!hasGwas && !hasTraits) return "";

    const { top_traits, ...gwasSlim } = gwasObj || {};
  
    const gwasPart = hasGwas ? `
      <details style="margin-top:6px;">
        <summary ${SECT_SUMMARY_STYLE}>GWAS profile</summary>
        <div style="margin-top:6px;">${objectToTable(gwasSlim, defaultEscapeHtml)}</div>
      </details>` : "";
  
    const traitsPart = hasTraits ? `
      <details style="margin-top:6px;">
        <summary ${SECT_SUMMARY_STYLE}>Traits (${traitsArr.length})</summary>
        <div style="margin-top:6px;">${traitsArr.map(t => `<span class="chip wrap">${escapeInline(String(t))}</span>`).join(" ")}</div>
      </details>` : "";
  
    return `
      <details style="margin-top:8px;">
        <summary ${TOP_SUMMARY_STYLE}>GWAS & Traits</summary>
        <div style="margin-top:6px;">
          ${gwasPart}
          ${traitsPart}
        </div>
      </details>
    `;
  }

  function renderInteractionsBlock({ sym, humanTSV, nonTSV, humanCount, nonhumanCount, goEnrichment }) {
    const inner = renderInteractionsInner({ sym, humanTSV, nonTSV, humanCount, nonhumanCount, goEnrichment });
    if (!inner) return "";
    return `
      <details style="margin-top:8px;">
        <summary ${TOP_SUMMARY_STYLE}>BioGRID Interactions</summary>
        <div style="margin-top:6px;">
          ${inner}
        </div>
      </details>
    `;
  }
  
  function renderInteractionsInner({ sym, humanTSV, nonTSV, humanCount, nonhumanCount, goEnrichment }) {
    const parts = [];
  
    // helper: parse TSV → { experiment_type -> Set(partner_symbol) }
    function groupInteractions(tsv) {
      const { header, rows } = parseTSVToTable(tsv);
      if (!header.length || !rows.length) return {};
      const idxType = header.findIndex(h => /experiment[_\s-]*type/i.test(h));
      const idxGene = header.findIndex(h => /(partner|interactor).*symbol/i.test(h));
      const map = {};
      rows.forEach(r => {
        const et = (r[idxType] || "").trim();
        const gene = (r[idxGene] || "").trim();
        if (!et || !gene) return;
        (map[et] ||= new Set()).add(gene);
      });
      return map;
    }
  
    function renderGrouped(kindLabel, tsv) {
      const groups = groupInteractions(tsv);
      const types = Object.keys(groups).sort((a,b) => a.localeCompare(b));
      if (!types.length) return "";
  
      const blocks = types.map(et => {
        const genes = Array.from(groups[et]).sort();
        const count = genes.length;
        const pills = renderGeneList(genes);
        return `
          <div style="margin:8px 0;">
            <div class="muted" style="font-weight:700; margin-bottom:4px;">${defaultEscapeHtml(et)} — ${count} gene${count===1?"":"s"}</div>
            <div class="row" style="gap:6px;flex-wrap:wrap;">${pills}</div>
          </div>`;
      }).join("");
  
      const totalPartners = types.reduce((a, et) => a + groups[et].size, 0);
      return `
        <details style="margin-top:8px;">
          <summary class="muted sect-hdr">
            ${defaultEscapeHtml(kindLabel)} (${types.length} types, ${totalPartners} partners)
          </summary>
          <div style="margin-top:6px;">${blocks}</div>
        </details>`;
    }
  
    // counts row (no extra "Interactions" header; parent summary already names section)
    if (humanCount > 0 || nonhumanCount > 0) {
      parts.push(`
        <div class="row" style="gap:6px; flex-wrap:wrap; margin-top:8px;">
          ${humanCount    > 0 ? `<span class="pill wrap">Human interactions: ${humanCount}</span>` : ""}
          ${nonhumanCount > 0 ? `<span class="pill wrap">Non-human interactions: ${nonhumanCount}</span>` : ""}
        </div>
      `);
    }
  
    if (humanCount > 0)    parts.push(renderGrouped("Human interactions", humanTSV));
    if (nonhumanCount > 0) parts.push(renderGrouped("Non-human interactions", nonTSV));
  
    // GO enrichment
    if (goEnrichment && (hasAny(goEnrichment.pan_go) || hasAny(goEnrichment.old_go))) {
      const blockId = `go_${sym}_${rand7()}`;
      const pan = goEnrichment.pan_go || {};
      const old = goEnrichment.old_go || {};
      parts.push(`
        <div class="muted" style="margin-top:8px;">GO Enrichment</div>
        <div class="row" style="gap:6px; margin-bottom:6px;">
          <label class="pill"><input type="radio" name="${blockId}" value="pan" checked> pan_go</label>
          <label class="pill"><input type="radio" name="${blockId}" value="old"> old_go</label>
        </div>
        <div id="${blockId}_content" data-go-block="${blockId}"
             data-go-pan='${jsonAttr(pan)}'
             data-go-old='${jsonAttr(old)}'
             class="mono" style="white-space:normal; line-height:1.35;">
          ${renderGoTriplet(pan, defaultEscapeHtml)}
        </div>
      `);
    }
  
    return parts.length ? parts.join("\n") : "";
  }
   
  

export async function renderGeneSummary(st, deps = {}) {
    const byId      = deps.byId || ((id) => document.getElementById(id));
    const show      = deps.show || ((id) => { const el = byId(id); if (el) el.style.display = "block"; });
    const hide      = deps.hide || ((id) => { const el = byId(id); if (el) el.style.display = "none"; });
    const fetchText = deps.fetchTextCached || (async (p) => (await fetch(p, { cache: "no-store" })).text());
    const parseCSV  = deps.parseCSV || defaultParseCSV;
    const esc       = deps.escapeHtml || defaultEscapeHtml;
  
    const mount = byId("geneSummary");
    const card  = byId("geneSummaryCard");
    if (!mount || !card) return;
  
    // collect symbols: STATE first, then index
    let idxRows = [];
    try { idxRows = parseCSV(await fetchText("genes/genes.index.csv"), ","); } catch {}
    const symKey  = guessKey(idxRows[0], ["symbol","gene_symbol","gene","SYMBOL","Gene"]);
    const fromState = Object.keys(st?.gene_entities || {}).map(safeUpper);
    const fromIndex = idxRows.map(r => safeUpper(r?.[symKey])).filter(Boolean);
    const symbols = uniquePreserve([...fromState, ...fromIndex]).filter(Boolean);
  
    if (!symbols.length) { mount.innerHTML = ""; hide("geneSummaryCard"); return; }
  
    const blocks = [];
    for (const sym of symbols) {
      const folder = `genes/${sym}`;
  
      // load gene.json (if missing, skip this gene)
      let g = null;
      try { g = JSON.parse(await fetchText(`${folder}/gene.json`) || "{}"); } catch {}
      if (!g || typeof g !== "object") continue;
  
      // optional TSVs for interaction counts/previews
      const humanTSV = await safeFetchText(fetchText, `${folder}/interactions_human.tsv`);
      const nonTSV   = await safeFetchText(fetchText, `${folder}/interactions_nonhuman.tsv`);
      const humanCount    = countDataRows(humanTSV);
      const nonhumanCount = countDataRows(nonTSV);
  
      // id header
      const id = g.identifiers || {};
      const idBits = [];
      if (id.ensembl_id) idBits.push(esc(id.ensembl_id));
      if (id.entrez_id)  idBits.push(`Entrez:${esc(id.entrez_id)}`);
      if (id.uniprot_id) idBits.push(`UniProt:${esc(id.uniprot_id)}`);
      const idLine = idBits.length ? ` • ${idBits.join(" • ")}` : "";
  
      // locus: use exactly as present (no inference)
      const loc = g.location || {};
      const locusHtml = hasLocus(loc) ? renderLocus(loc, esc) : "";
  
      // transcript stats: only if transcript_count > 0 (or other non-null stats)
      const tx = g.transcriptome || {};
      const txBits = [];
      if (isPosInt(tx.transcript_count)) txBits.push(`Transcripts: ${esc(String(tx.transcript_count))}`);
      if (isNum(tx.median_transcript_span_bp)) txBits.push(`Median span: ${esc(String(tx.median_transcript_span_bp))} bp`);
      if (isNum(tx.max_transcript_span_bp))    txBits.push(`Max span: ${esc(String(tx.max_transcript_span_bp))} bp`);
      const txLine = txBits.length ? txBits.join(" • ") : "";
  
      // expandable sections (only when non-empty)
      const sections = [];
  
      // Annotations (arrays → pills/table-ish)
      sections.push(renderArraySection("Functions", g.annotations?.functions));
      sections.push(renderArraySection("Diseases",  g.annotations?.diseases));
      sections.push(renderArraySection("Pathways",  g.annotations?.pathways));
      sections.push(renderArraySection("Cell localization", g.annotations?.cell_localizations));
      sections.push(renderArraySection("GO annotations",   g.annotations?.go_annotations));
      sections.push(renderArraySection("Predicted GO",     g.annotations?.predicted_go));
      sections.push(renderArraySection("Predicted disease",g.annotations?.predicted_disease));
  
      // MSigDB + OMIM
      sections.push(renderMsigdbBlock(g.msigDB?.msigdb_annotations));
// OMIM 
sections.push(renderArraySection("OMIM phenotypes", g.omim?.phenotype_annotations));

// OpenTargets (disease dict/array + other fields + conditional is_essential)
sections.push(renderOpenTargetsBlock(g.open_targets));
      sections.push(renderObjectSection("Tissue-specific expression", g.open_targets?.tissue_specific_expression));
  
      // Genetic constraint (object)
      sections.push(renderObjectSection("Genetic constraint", g.genetic_constraint));
  
      // Transcriptome sub-objects (if they’re not empty)
      sections.push(renderTranscriptsSection(tx.transcripts));
sections.push(renderIsoformsSection(tx.isoforms));

  
      // GWAS profile (object) — if present
      sections.push(renderGwasTraitsBlock(g.gwas_profile, g.traits));
      
      sections.push(renderObjectSection("Binding peaks", g.binding_peaks));
      sections.push(renderObjectSection("miRNA targets", g.mirna_targets));
  
      // Interactions: TSV counts + collapsible previews + special GO enrichment toggle
      const interactionsBlock = renderInteractionsBlock({
        sym, humanTSV, nonTSV, humanCount, nonhumanCount, goEnrichment: g.interactions?.go_enrichment
      });
  
      // assemble card
      blocks.push(`
        <div class="card" style="margin-bottom:10px;">
          <div style="display:flex;align-items:center;gap:8px;flex-wrap:wrap;">
            <strong style="font-size:18px;">${esc(sym)}</strong>
            <span class="muted">${idLine}</span>
          </div>
  
          ${ (locusHtml || txLine) ? `
            <div class="muted" style="margin-top:4px;">
              ${locusHtml}${locusHtml && txLine ? " • " : ""}${txLine}
            </div>` : "" }
  
          ${sections.filter(Boolean).join("\n")}
          ${interactionsBlock}
        </div>
      `);
    }
  
    if (!blocks.length) { mount.innerHTML = ""; hide("geneSummaryCard"); return; }
    mount.innerHTML = blocks.join("\n");
    wireExpanders(mount);
    show("geneSummaryCard");
  
    /* ----- helpers (render) ----- */

    function isGeneralLocalization(o){
        // matches both: object itself is a general_localization, or nested under { general_localization: {…} }
        if (!o || typeof o !== "object") return false;
        if (o.name === "general_localization" || o.general_localization) return true;
        return false;
      }
      
      function coerceGeneralLocalization(o){
        return o.name === "general_localization" ? o : o.general_localization || null;
      }
      
      function renderGeneralLocalizationSection(title, obj){
        const gl = coerceGeneralLocalization(obj);
        if (!gl) return "";
        const locs  = Array.isArray(gl.locations) ? gl.locations : [];
        const notes = Array.isArray(gl.notes) ? gl.notes : [];
      
        const locRow  = locs.length  ? `<div class="row" style="gap:6px;flex-wrap:wrap;">${locs.map(x=>`<span class="chip">${esc(String(x))}</span>`).join("")}</div>` : "";
        // long notes → tiny chips that wrap; big blocks are noisy
        const noteRow = notes.length ? `<div class="row" style="gap:6px;flex-wrap:wrap;margin-top:6px;">${notes.map(n=>`<span class="chip chip-note">${esc(String(n))}</span>`).join("")}</div>` : "";
      
        const meta = [];
        if (gl.status) meta.push(`Status: <strong>${esc(String(gl.status))}</strong>`);
        if (Array.isArray(gl.aliases) && gl.aliases.length) meta.push(`Aliases: ${esc(gl.aliases.join(", "))}`);
      
        return `
          <details style="margin-top:10px;">
            <summary class="muted sect-hdr">${esc(title)} — General localization</summary>
            <div style="margin-top:8px;">
              ${meta.length ? `<div class="muted" style="margin-bottom:6px;">${meta.join(" • ")}</div>` : ""}
              ${locRow || `<div class="muted">No locations reported.</div>`}
              ${noteRow}
            </div>
          </details>
        `;
      }
      

      function renderArraySection(title, arr) {
        if (!Array.isArray(arr) || arr.length === 0) return "";
      
        // SPECIAL CASE: Functions -> plain wrapped text (no chips)
        if (title === "Functions") {
          // Join items into a single readable block. 
          const text = arr.map(x => String(x).trim()).filter(Boolean).join(" ");
          return `
            <details style="margin-top:8px;">
              <summary ${TOP_SUMMARY_STYLE}>
                ${escapeInline(title)} (${arr.length})
              </summary>
              <div class="blk-text" style="margin-top:6px; white-space:normal; line-height:1.4; word-break:break-word; overflow-wrap:anywhere;">
                ${escapeInline(text)}
              </div>
            </details>
          `;
        }
      
        // Default: render as chips
        const list = arr.map(x => `<span class="chip wrap">${escapeInline(String(x))}</span>`).join(" ");
        return `
          <details style="margin-top:6px;">
            <summary ${TOP_SUMMARY_STYLE}>
              ${escapeInline(title)} (${arr.length})
            </summary>
<div class="row" style="margin-top:4px;">${list || `<div class="muted">No items.</div>`}</div>
          </details>
        `;
      }
      
      
      function renderObjectSection(title, obj) {
        if (!obj || typeof obj !== "object" || Object.keys(obj).length === 0) return "";
        const table = objectToTable(obj, esc);
        return `
          <details style="margin-top:6px;">
            <summary ${TOP_SUMMARY_STYLE}>
              ${esc(title)}
            </summary>
            <div style="margin-top:6px;">${table}</div>
          </details>
        `;
      }
      }
      
  
  /* ---------------- behaviors ---------------- */
  
  function wireExpanders(root) {
    // Lazy TSV mini-tables when opening details
    root.querySelectorAll('details[data-lazy-table]').forEach(d => {
      let rendered = false;
      d.addEventListener('toggle', () => {
        if (d.open && !rendered) {
          rendered = true;
          try {
            const cfg = JSON.parse(d.getAttribute('data-lazy-table') || "{}");
            const target = root.querySelector(`#${cfg.id}`);
            if (target) {
              const { header, rows } = parseTSVToTable(cfg.tsv);
              target.innerHTML = renderMiniTable(header, rows.slice(0, cfg.max || 30));
              if (rows.length > (cfg.max || 30)) {
                const more = document.createElement('div');
                more.className = 'muted';
                more.style.marginTop = '6px';
                more.textContent = `Showing first ${cfg.max || 30} of ${rows.length} rows.`;
                target.appendChild(more);
              }
            }
          } catch {}
        }
      });
    });
  
    // GO enrichment toggles
    root.querySelectorAll('[data-go-block]').forEach(el => {
      const blockId = el.getAttribute('data-go-block');
      const radios = root.querySelectorAll(`input[name="${blockId}"]`);
      const pan = tryJSON(el.getAttribute('data-go-pan')) || {};
      const old = tryJSON(el.getAttribute('data-go-old')) || {};
      const render = (obj) => { el.innerHTML = renderGoTriplet(obj, defaultEscapeHtml); };
      radios.forEach(r => r.addEventListener('change', () => render(r.value === "pan" ? pan : old)));
    });
  }
  
  /* ---------------- tiny render helpers ---------------- */
  
  function objectToTable(obj, esc){
    // flatten simple objects; arrays render as comma-joined; nested objects -> JSON snippet
    const rows = Object.entries(obj).map(([k,v]) => {
      let val = "";
      if (Array.isArray(v)) val = v.map(x => String(x)).join(", ");
      else if (v && typeof v === "object") val = JSON.stringify(v);
      else val = String(v ?? "");
      return `<tr><th class="wrap" style="text-align:left;vertical-align:top;">${esc(k)}</th><td class="wrap">${esc(val)}</td></tr>`;
    }).join("");
    return `<div style="overflow:auto; max-height:300px;"><table class="vs-table"><tbody>${rows}</tbody></table></div>`;
  }
  
  function renderMiniTable(header, rows){
    const thead = `<thead><tr>${header.map(h => `<th class="wrap">${escapeInline(h)}</th>`).join("")}</tr></thead>`;
    const tbody = `<tbody>${rows.map(r => `<tr>${r.map(v => `<td class="wrap">${escapeInline(v)}</td>`).join("")}</tr>`).join("")}</tbody>`;
    return `<div style="overflow:auto; max-height:300px;"><table class="vs-table" style="table-layout:auto;width:100%;">${thead}${tbody}</table></div>`;
  }
  
  function renderGoTriplet(obj, esc){
    const sect = ["BP","MF","CC"].map(k => {
      const arr = Array.isArray(obj[k]) ? obj[k] : [];
      if (!arr.length) return "";
      return `<div style="margin-top:4px;"><strong>${k}</strong>: ${arr.slice(0,8).map(x => esc(String(x))).join("; ")}</div>`;
    }).join("");
    return sect || `<div class="muted">No GO enrichment.</div>`;
  }
  
  /* ---------------- parsing / checks ---------------- */
  
  function defaultEscapeHtml(s){ return String(s ?? "").replace(/[&<>"']/g, m => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[m])); }
  function escapeInline(s){ return defaultEscapeHtml(String(s ?? "")); }
  
  function defaultParseCSV(text, delimiter = ","){
    if (!text) return [];
    const lines = String(text).replace(/\r\n?/g, "\n").split("\n").filter(l => l.trim());
    if (!lines.length) return [];
    const parseLine = (line) => {
      const out = []; let cur = "", inQ = false;
      for (let i=0;i<line.length;i++){
        const ch = line[i];
        if (inQ){
          if (ch === '"' && line[i+1] === '"'){ cur+='"'; i++; }
          else if (ch === '"'){ inQ = false; }
          else cur += ch;
        } else {
          if (ch === '"') inQ = true;
          else if (ch === delimiter){ out.push(cur); cur=""; }
          else cur += ch;
        }
      }
      out.push(cur);
      return out;
    };
    const header = parseLine(lines[0]).map(h => h.trim());
    return lines.slice(1).map(raw => {
      const parts = parseLine(raw);
      const obj = {};
      for (let i=0;i<header.length;i++) obj[header[i]] = (parts[i] ?? "").trim();
      return obj;
    });
  }
  
  function parseTSVToTable(tsv){
    if (!tsv) return { header: [], rows: [] };
    const lines = String(tsv).split(/\r?\n/).filter(l => l.trim() && !/^\s*#/.test(l));
    if (!lines.length) return { header: [], rows: [] };
    const header = lines[0].split("\t").map(s => s.trim());
    const rows = lines.slice(1).map(l => l.split("\t"));
    const hasData = rows.some(r => r.some(x => String(x||"").trim()));
    return hasData ? { header, rows } : { header, rows: [] };
  }
  
  function countDataRows(tsv){
    const { rows } = parseTSVToTable(tsv);
    return rows.length;
  }
  
  function tryJSON(s){ try { return JSON.parse(s || ""); } catch { return null; } }
  
  /* ---------------- misc utils ---------------- */
  
  function hasLocus(loc){
    if (!loc || typeof loc !== "object") return false;
    return !!(loc.chrom || loc.chromosome || loc.start != null || loc.end != null || loc.strand);
  }
  function renderLocus(loc, esc){
    const chrom = String(loc.chrom || loc.chromosome || "").trim();
    const s = loc.start != null ? String(loc.start) : "?";
    const e = loc.end   != null ? String(loc.end)   : "?";
    const st = loc.strand ? ` ${esc(String(loc.strand))}` : "";
    if (!chrom && s==="?" && e==="?") return "";
    const c = chrom ? `chr${esc(stripChr(chrom))}` : "chr?";
    return `${c} ${s}–${e}${st}`;
  }
  function stripChr(x){ return String(x).replace(/^chr/i,"").replace(/^m(t)?$/i,"MT"); }
  
  function hasAny(x){
    if (!x) return false;
    if (Array.isArray(x)) return x.length > 0;
    if (typeof x === "object") return Object.keys(x).length > 0;
    return !!x;
  }
  function isPosInt(x){ return Number.isInteger(x) && x > 0; }
  function isNum(x){ const n = Number(x); return Number.isFinite(n); }
  
  function safeUpper(s){ s = (s==null?"":String(s)).trim(); return s ? s.toUpperCase() : ""; }
  function uniquePreserve(a){ const seen=new Set(); const out=[]; for(const x of a){ if(!seen.has(x)){ seen.add(x); out.push(x);} } return out; }
  function guessKey(row, cands){ if(!row||typeof row!=="object") return null; const keys=Object.keys(row); for(const w of cands){ const ex=keys.find(k=>k.toLowerCase()===String(w).toLowerCase()); if(ex) return ex; } const lk=keys.map(k=>[k,k.toLowerCase()]); for(const w of cands){ const s=String(w).toLowerCase(); const hit=lk.find(([k,l])=>l.includes(s)); if(hit) return hit[0]; } return null; }
  function rand7(){ return Math.random().toString(36).slice(2,9); }
  
  async function safeFetchText(fetcher, path){ try { return await fetcher(path); } catch { return ""; } }
  function jsonAttr(o){ try { return defaultEscapeHtml(JSON.stringify(o)); } catch { return "{}"; } }
  