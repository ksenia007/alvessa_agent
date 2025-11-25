// src/tools/aa_seq/frontend/tool_aa_seq.js
// Author: Dmitri Kosenkov
// Updated: 2025-11-25
//
// AA sequence -> UniProt / gene mapping viewer.
//
// Expected global injected by Python (from aa_seq.utils.inject_frontend_assets_aa_seq):
// const aaSeqData = {
//   genes: {
//     CRIPTO: {
//       entrez_gene_id: "6997",
//       n_hits: 1,
//       hits: [
//         {
//           acc: "P13385",
//           canonical_acc: "P13385",
//           score: 100.0,
//           coverage: 1.0,
//           alignment_length: 30,
//           query_sequence: "CKCWHGQLRCFPQAFLPGCDGLVMDEHLVA",
//           uniprot_sequence: "MDCRKMARF...LSIQSYY",
//           source: ""
//         },
//         ...
//       ]
//     },
//     ...
//   },
//   unmapped: [ { ... same shape, no gene_name ... } ],
//   meta: {
//     total_hits: 2,
//     distinct_genes: 2,
//     total_unmapped: 0,
//     sequences: ["CKCWHGQLRCFPQAFLPGCDGLVMDEHLVA", ...]
//   }
// };

let currentGene = null;

// =============================
// Helpers
// =============================
function esc(s) {
  if (s == null) return "";
  return String(s)
    .replace(/&/g, "&amp;")
    .replace(/</g, "&lt;")
    .replace(/>/g, "&gt;")
    .replace(/"/g, "&quot;")
    .replace(/'/g, "&#39;");
}

function formatScore(value) {
  if (value == null || value === "" || Number.isNaN(Number(value))) {
    return "N/A";
  }
  const v = Number(value);
  return v.toFixed(2);
}

function formatCoverage(value) {
  if (value == null || value === "" || Number.isNaN(Number(value))) {
    return "N/A";
  }
  const v = Number(value);
  return v.toFixed(3);
}

function formatAlignLen(value) {
  if (value == null || value === "" || Number.isNaN(Number(value))) {
    return "N/A";
  }
  return String(value);
}

function truncateSeq(seq, maxLen) {
  if (!seq) return "";
  const s = String(seq);
  if (s.length <= maxLen) return s;
  return s.slice(0, maxLen) + "...";
}

function getQuerySeq(hit) {
  // support both "sequence" and "query_sequence" just in case
  if (!hit) return "";
  const v = hit.sequence || hit.query_sequence || "";
  return String(v).trim();
}

function isRealGeneKey(key) {
  if (!key) return false;
  // we no longer use __meta__/__unmatched__ keys, but keep guard
  return !key.startsWith("__");
}

// =============================
// Table builders
// =============================
function buildMappedTable(tbody, geneName, geneData) {
  if (!tbody) return;
  tbody.innerHTML = "";

  if (!geneData || !Array.isArray(geneData.hits) || geneData.hits.length === 0) {
    return;
  }

  const hits = geneData.hits;
  if (!hits.length) return;

  for (const hit of hits) {
    const tr = document.createElement("tr");

    const q = getQuerySeq(hit);
    const ref = (hit.uniprot_sequence || "").trim();
    const acc = hit.acc || "";
    const canon = hit.canonical_acc || "";
    const egid = geneData.entrez_gene_id ? String(geneData.entrez_gene_id) : "";

    // Gene / Accession
    const tdGeneAcc = document.createElement("td");
    let inner = esc(geneName || "-");
    if (acc) {
      inner += '<div class="muted">' + esc(acc) + "</div>";
    }
    if (canon && canon !== acc) {
      inner += '<div class="muted">Canonical: ' + esc(canon) + "</div>";
    } else if (canon && canon === acc) {
      inner += '<div class="muted">Canonical accession</div>';
    }
    tdGeneAcc.innerHTML = inner;
    tr.appendChild(tdGeneAcc);

    // Entrez Gene ID
    const tdEntrez = document.createElement("td");
    tdEntrez.textContent = egid || "-";
    tr.appendChild(tdEntrez);

    // Score (%)
    const tdScore = document.createElement("td");
    tdScore.className = "num";
    tdScore.textContent = formatScore(hit.score);
    tr.appendChild(tdScore);

    // Coverage
    const tdCov = document.createElement("td");
    tdCov.className = "num";
    tdCov.textContent = formatCoverage(hit.coverage);
    tr.appendChild(tdCov);

    // Alignment length
    const tdLen = document.createElement("td");
    tdLen.className = "num";
    tdLen.textContent = formatAlignLen(hit.alignment_length);
    tr.appendChild(tdLen);

    // Query length
    const tdQLen = document.createElement("td");
    tdQLen.className = "num";
    tdQLen.textContent = q ? String(q.length) : "N/A";
    tr.appendChild(tdQLen);

    // Query sequence
    const tdQuery = document.createElement("td");
    tdQuery.className = "seq-cell";
    tdQuery.textContent = q;
    tdQuery.title = q;
    tr.appendChild(tdQuery);

    // UniProt reference sequence
    const tdRef = document.createElement("td");
    tdRef.className = "seq-cell";
    tdRef.textContent = ref;
    tdRef.title = ref;
    tr.appendChild(tdRef);

    tbody.appendChild(tr);
  }
}

function buildUnmappedTable(tbody, unmatched) {
  if (!tbody) return;
  tbody.innerHTML = "";

  if (!Array.isArray(unmatched) || unmatched.length === 0) return;

  for (const hit of unmatched) {
    const tr = document.createElement("tr");

    const q = getQuerySeq(hit);
    const ref = (hit.uniprot_sequence || "").trim();
    const acc = hit.acc || "";
    const canon = hit.canonical_acc || "";

    // Accession
    const tdAcc = document.createElement("td");
    tdAcc.textContent = acc || "-";
    tr.appendChild(tdAcc);

    // Canonical accession
    const tdCanon = document.createElement("td");
    if (canon && canon !== acc) {
      tdCanon.textContent = canon;
    } else if (canon && canon === acc) {
      tdCanon.textContent = canon;
    } else {
      tdCanon.textContent = "-";
    }
    tr.appendChild(tdCanon);

    // Score (%)
    const tdScore = document.createElement("td");
    tdScore.className = "num";
    tdScore.textContent = formatScore(hit.score);
    tr.appendChild(tdScore);

    // Coverage
    const tdCov = document.createElement("td");
    tdCov.className = "num";
    tdCov.textContent = formatCoverage(hit.coverage);
    tr.appendChild(tdCov);

    // Alignment length
    const tdLen = document.createElement("td");
    tdLen.className = "num";
    tdLen.textContent = formatAlignLen(hit.alignment_length);
    tr.appendChild(tdLen);

    // Query length
    const tdQLen = document.createElement("td");
    tdQLen.className = "num";
    tdQLen.textContent = q ? String(q.length) : "N/A";
    tr.appendChild(tdQLen);

    // Query sequence
    const tdQuery = document.createElement("td");
    tdQuery.className = "seq-cell";
    tdQuery.textContent = q;
    tdQuery.title = q;
    tr.appendChild(tdQuery);

    // UniProt reference sequence
    const tdRef = document.createElement("td");
    tdRef.className = "seq-cell";
    tdRef.textContent = ref;
    tdRef.title = ref;
    tr.appendChild(tdRef);

    tbody.appendChild(tr);
  }
}

// =============================
// Populate view
// =============================
function populateGene(geneName) {
  if (typeof aaSeqData === "undefined" || !aaSeqData) {
    return;
  }

  currentGene = geneName || null;

  const mappedSection = document.getElementById("mappedSection");
  const unmappedSection = document.getElementById("unmappedSection");

  if (!mappedSection || !unmappedSection) {
    // Template mismatch; nothing we can safely do
    return;
  }

  const mappedHeader = mappedSection.querySelector(".collapsible-header .title");
  const mappedCount = mappedSection.querySelector(".collapsible-header .count");
  const mappedBody = document.getElementById("mappedTable");

  const unmappedHeader = unmappedSection.querySelector(".collapsible-header .title");
  const unmappedCount = unmappedSection.querySelector(".collapsible-header .count");
  const unmappedBody = document.getElementById("unmappedTable");

  const meta = (aaSeqData && aaSeqData.meta) || {};
  const genesDict = (aaSeqData && aaSeqData.genes) || {};
  const unmatched = (aaSeqData && aaSeqData.unmapped) || [];

  // Per-gene hits
  const geneData = geneName ? genesDict[geneName] : null;
  if (!geneData || !Array.isArray(geneData.hits) || geneData.hits.length === 0) {
    mappedSection.style.display = "none";
  } else {
    mappedSection.style.display = "block";
    const nHits = geneData.hits.length;

    const egid = geneData.entrez_gene_id ? String(geneData.entrez_gene_id) : "";
    if (mappedHeader) {
      if (egid) {
        mappedHeader.textContent =
          "Mapped UniProt hits for " + geneName + " (Entrez Gene ID: " + egid + ")";
      } else {
        mappedHeader.textContent = "Mapped UniProt hits for " + geneName;
      }
    }
    if (mappedCount) mappedCount.textContent = nHits;

    buildMappedTable(mappedBody, geneName, geneData);
  }

  // Unmapped section (accessions with no gene_name)
  if (!unmatched.length) {
    unmappedSection.style.display = "none";
  } else {
    unmappedSection.style.display = "block";
    const nUnmatched = unmatched.length;

    if (unmappedHeader) {
      unmappedHeader.textContent =
        "UniProt accessions without gene annotation (all sequences)";
    }
    if (unmappedCount) unmappedCount.textContent = nUnmatched;

    buildUnmappedTable(unmappedBody, unmatched);
  }

  // Optional debug
  if (
    meta &&
    typeof meta.total_hits === "number" &&
    typeof meta.distinct_genes === "number"
  ) {
    // eslint-disable-next-line no-console
    console.log(
      "[AA-SEQ] total hits = " +
        meta.total_hits +
        ", distinct genes = " +
        meta.distinct_genes +
        ", total_unmapped = " +
        (meta.total_unmapped != null ? meta.total_unmapped : "N/A")
    );
  }
}

// =============================
// Initialize selects
// =============================
function initGeneSelect() {
  if (typeof aaSeqData === "undefined" || !aaSeqData) return;

  const geneSelect = document.getElementById("geneSelect");
  if (!geneSelect) return;

  const genesDict = aaSeqData.genes || {};
  const genes = [];

  Object.keys(genesDict).forEach((key) => {
    if (!isRealGeneKey(key)) return;
    const gd = genesDict[key];
    const hits = gd && Array.isArray(gd.hits) ? gd.hits.length : 0;
    if (hits > 0) {
      genes.push({ name: key, hits });
    }
  });

  genes.sort((a, b) => a.name.localeCompare(b.name));

  if (!genes.length) {
    geneSelect.style.display = "none";
    const mappedSection = document.getElementById("mappedSection");
    if (mappedSection) mappedSection.style.display = "none";
    return;
  }

  for (const g of genes) {
    const opt = document.createElement("option");
    opt.value = g.name;
    opt.textContent = g.name + " (" + g.hits + ")";
    geneSelect.appendChild(opt);
  }

  currentGene = genes[0].name;
  geneSelect.value = currentGene;
  geneSelect.addEventListener("change", (e) => {
    currentGene = e.target.value || null;
    populateGene(currentGene);
  });
}

// =============================
// Collapsibles + theme sync
// =============================
function initCollapsibles() {
  document.querySelectorAll(".collapsible-header").forEach((header) => {
    const section = header.closest(".collapsible-section");
    if (!section) return;
    const content = section.querySelector(".collapsible-content");
    const icon = header.querySelector(".toggle-icon");
    if (!content || !icon) return;

    header.addEventListener("click", () => {
      const isOpen = content.classList.toggle("open");
      header.classList.toggle("active", isOpen);
      icon.textContent = isOpen ? "▾" : "▸";
      content.style.maxHeight = "";
    });

    // start collapsed
    content.classList.remove("open");
    content.style.maxHeight = "";
    icon.textContent = "▸";
  });
}

// =============================
// DOM Ready
// =============================
document.addEventListener("DOMContentLoaded", () => {
  initGeneSelect();
  initCollapsibles();

  // Initial populate (if any genes exist)
  if (currentGene) {
    populateGene(currentGene);
  } else {
    // at least populate unmatched table if present
    populateGene(null);
  }

  // Secure theme sync (no same-origin required)
  window.addEventListener("message", (e) => {
    if (!e.data || typeof e.data.theme !== "string") return;
    if (e.data.theme === "light") {
      document.body.classList.add("light");
    } else {
      document.body.classList.remove("light");
    }
  });
});
