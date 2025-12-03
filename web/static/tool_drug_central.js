// web/static/tool_drug_central.js
// Drug Central drug-centric viewer
// Expects a single global:
//   var drug_centralData = { "DC:3290": { ... }, ... };

let currentDrug = null;
let RDKit = null;

// =============================
// Helpers
// =============================
function drug_centralUrl(structId) {
  if (structId == null || structId === "") return "https://drugcentral.org/";
  return "https://drugcentral.org/drugcard/" + encodeURIComponent(String(structId));
}

function pubchemCidUrl(cid) {
  return "https://pubchem.ncbi.nlm.nih.gov/compound/" + cid;
}

function esc(s) {
  if (s == null) return "";
  return String(s)
    .replace(/&/g, "&amp;")
    .replace(/</g, "&lt;")
    .replace(/>/g, "&gt;")
    .replace(/"/g, "&quot;")
    .replace(/'/g, "&#39;");
}

// Map Drug Central id_type -> nice label
function mapIdType(idType) {
  if (!idType) return "Other";
  switch (idType.toUpperCase()) {
    case "DRUG_CENTRAL_ID":
      return "DrugCentral";
    case "STRUCT_ID":
      return "Struct ID";
    case "PUBCHEM_CID":
      return "PubChem CID";
    case "UNII":
      return "UNII";
    case "CHEBI":
      return "ChEBI";
    case "CHEMBL_ID":
    case "CHEMBL":
      return "ChEMBL ID";
    case "RXNORM":
      return "RxNorm";
    default:
      return idType;
  }
}

// =============================
// PubChem Lookup (cached + throttled)
// =============================
const _pubchemCache = new Map(); // inchikey -> { cid, title }

async function fetchPubChemInfo(inchikey) {
  if (!inchikey) return null;
  if (_pubchemCache.has(inchikey)) return _pubchemCache.get(inchikey);
  try {
    // simple throttle
    await new Promise(function (res) { setTimeout(res, 200); });

    const cidResp = await fetch(
      "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/" +
        encodeURIComponent(inchikey) +
        "/cids/JSON"
    );
    if (!cidResp.ok) {
      console.warn("PubChem CID status " + cidResp.status + " for " + inchikey);
      _pubchemCache.set(inchikey, null);
      return null;
    }
    const cidData = await cidResp.json();
    if (!cidData.IdentifierList || !cidData.IdentifierList.CID || !cidData.IdentifierList.CID.length) {
      _pubchemCache.set(inchikey, null);
      return null;
    }
    const cid = cidData.IdentifierList.CID[0];

    await new Promise(function (res) { setTimeout(res, 200); });
    const propResp = await fetch(
      "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/" +
        cid +
        "/property/Title/JSON"
    );
    if (!propResp.ok) {
      console.warn("PubChem property status " + propResp.status + " for " + inchikey);
      _pubchemCache.set(inchikey, null);
      return null;
    }
    const propData = await propResp.json();
    const title =
      propData &&
      propData.PropertyTable &&
      propData.PropertyTable.Properties &&
      propData.PropertyTable.Properties[0] &&
      propData.PropertyTable.Properties[0].Title
        ? propData.PropertyTable.Properties[0].Title
        : "";

    const out = { cid: cid, title: title };
    _pubchemCache.set(inchikey, out);
    return out;
  } catch (err) {
    console.warn("PubChem lookup error:", err, "for", inchikey);
    _pubchemCache.set(inchikey, null);
    return null;
  }
}

// =============================
// RDKit SVG rendering
// =============================
function renderMol2D(cell, smiles, inchikey) {
  if (!RDKit || !smiles) {
    cell.textContent = "-";
    return;
  }
  try {
    const mol = RDKit.get_mol(smiles);
    if (!mol) throw new Error("RDKit mol creation failed");
    const svg = mol.get_svg();
    mol.delete();

    const wrapper = document.createElement("div");
    wrapper.className = "mol-thumb";
    wrapper.innerHTML = svg;

    const svgElem = wrapper.querySelector("svg");
    if (svgElem) {
      svgElem.style.maxWidth = "90px";
      svgElem.style.cursor = "zoom-in";
    }

    wrapper.addEventListener("click", function () {
      openMolModal(svg, inchikey, smiles);
    });

    cell.innerHTML = "";
    cell.appendChild(wrapper);
  } catch (err) {
    console.warn("RDKit render failed:", err);
    cell.textContent = "-";
  }
}

// =============================
// Modal Viewer
// =============================
async function openMolModal(svg, inchikey, smiles) {
  const modal = document.createElement("div");
  modal.className = "mol-modal";
  modal.innerHTML = [
    '<div class="mol-modal-backdrop"></div>',
    '<div class="mol-modal-content">',
    '  <div class="mol-modal-header">',
    '    <span>Molecule View</span>',
    '    <button class="mol-modal-close" aria-label="Close">X</button>',
    "  </div>",
    '  <div class="mol-modal-body">' + svg + "</div>",
    '  <div class="mol-modal-footer">',
    '    <span id="molName">Loading PubChem...</span>',
    '    <button id="saveSvgBtn">Save SVG</button>',
    "  </div>",
    "</div>"
  ].join("");

  document.body.appendChild(modal);
  requestAnimationFrame(function () {
    modal.classList.add("active");
  });

  modal.querySelector(".mol-modal-close").onclick = function () { modal.remove(); };
  modal.querySelector(".mol-modal-backdrop").onclick = function () { modal.remove(); };
  modal.querySelector("#saveSvgBtn").onclick = function () {
    const blob = new Blob([svg], { type: "image/svg+xml" });
    const a = document.createElement("a");
    a.href = URL.createObjectURL(blob);
    a.download = (inchikey || "molecule") + ".svg";
    a.click();
  };

  const info = await fetchPubChemInfo(inchikey);
  const molNameElem = modal.querySelector("#molName");
  if (info && molNameElem) {
    const cid = info.cid;
    const title = info.title;
    if (title) {
      molNameElem.innerHTML =
        '<a href="' +
        pubchemCidUrl(cid) +
        '" target="_blank" rel="noopener">' +
        esc(title) +
        "</a>";
    } else if (cid) {
      molNameElem.innerHTML =
        '<a href="' +
        pubchemCidUrl(cid) +
        '" target="_blank" rel="noopener">PubChem CID ' +
        cid +
        "</a>";
    } else {
      molNameElem.style.display = "none";
    }
  } else if (molNameElem) {
    molNameElem.style.display = "none";
  }
}

// =============================
// Section builders
// =============================
function buildOverviewCard(container, data, drugKey) {
  container.innerHTML = "";
  if (!data) {
    container.textContent = "No overview data available.";
    return;
  }

  const dg = data.drug || {};
  const cr = data.cross_refs || {};

  const primaryName = dg.name || drugKey || "";
  const structId = dg.struct_id != null ? dg.struct_id : "";

  const smiles = dg.smiles || "";
  const inchikey = dg.inchikey || "";
  const cas = dg.cas_number || "";
  const status = dg.status || "";
  const atcCodes = (cr.atc_codes || []).filter(Boolean);
  const synonyms = (dg.synonyms || []).filter(Boolean);

  const left = document.createElement("div");
  left.className = "drug_central-card-left";

  const right = document.createElement("div");
  right.className = "drug_central-card-right";

  const dl = document.createElement("dl");

  function addField(label, valueHtml) {
    if (!valueHtml) return;
    const dt = document.createElement("dt");
    dt.textContent = label;
    const dd = document.createElement("dd");
    dd.innerHTML = valueHtml;
    dl.appendChild(dt);
    dl.appendChild(dd);
  }

  addField("Primary name", esc(primaryName));

  // Use DrugCentral struct_id (structures.id) for the external link
  if (structId !== "") {
    addField(
      "DrugCentral",
      '<a href="' +
        drug_centralUrl(structId) +
        '" target="_blank" rel="noopener">' +
        esc(String(structId)) +
        "</a>"
    );
  }

  //if (structId !== "") addField("Struct ID", esc(structId));
  if (cas) addField("CAS Number", esc(cas));
  if (status) addField("Status", esc(status));
  if (atcCodes && atcCodes.length) {
    addField("ATC codes", esc(atcCodes.join(", ")));
  }
  if (synonyms && synonyms.length) {
    addField("Synonyms", esc(synonyms.slice(0, 8).join(", ")));
  }
  if (inchikey) addField("InChIKey", esc(inchikey));
  if (smiles) addField("SMILES", "<code>" + esc(smiles) + "</code>");

  left.appendChild(dl);

  const molContainer = document.createElement("div");
  right.appendChild(molContainer);
  if (smiles) {
    renderMol2D(molContainer, smiles, inchikey);
  } else {
    molContainer.textContent = "No structure available.";
  }

  const wrapper = document.createElement("div");
  wrapper.className = "drug_central-card";
  wrapper.appendChild(left);
  wrapper.appendChild(right);
  container.appendChild(wrapper);
}

function buildRegulatoryTables(regBody, prodBody, data) {
  regBody.innerHTML = "";
  prodBody.innerHTML = "";
  if (!data) return { approvalsCount: 0, productsCount: 0 };

  const reg = data.regulatory || {};
  const approvals = reg.approvals || [];
  const obProducts = reg.orange_book_products || [];
  const regStatus = reg.status || "";

  // Approvals -> first table
  // Expected HTML header:
  //   Approval Date | Agency / Type | Applicant | Orphan | Regulatory Status | Notes
  approvals.forEach(function (row) {
    const tr = document.createElement("tr");
    function td(text) {
      const cell = document.createElement("td");
      cell.textContent = text || "-";
      return cell;
    }

    const dateStr =
      row.approval ||      // approval.approval (text / date string)
      row.approval_date || // generic date
      row.date ||          // fallback
      "";

    const agencyType =
      row.agency ||        // explicit agency if Python added it
      row.type ||          // approval.type (e.g. FDA, EMA)
      "Unknown";

    const applicant =
      row.applicant ||     // approval.applicant
      row.company ||       // alt naming
      "";

    const orphanFlag =
      row.orphan === 1 ||
      row.orphan === true ||
      row.orphan === "Y"
        ? "Yes"
        : "";

    const notesPieces = [];
    if (row.source) notesPieces.push(String(row.source));
    if (row.appl_type || row.appl_no) {
      notesPieces.push(
        [row.appl_type, row.appl_no].filter(Boolean).join(" ")
      );
    }
    if (row.notes) notesPieces.push(String(row.notes));
    const notes = notesPieces.join(" | ");

    tr.appendChild(td(dateStr));
    tr.appendChild(td(agencyType));
    tr.appendChild(td(applicant));
    tr.appendChild(td(orphanFlag));
    tr.appendChild(td(regStatus || ""));
    tr.appendChild(td(notes));
    regBody.appendChild(tr);
  });

  // Orange Book products -> second table
  // Expected HTML header:
  //   Product (Trade Name) | Ingredient(s) | Route | Dose Form | Strength |
  //   Application (Type / No. / Prod.) | Approval Date | Product Type / TE Code
  obProducts.forEach(function (row) {
    const tr = document.createElement("tr");
    function td(text) {
      const cell = document.createElement("td");
      cell.textContent = text || "-";
      return cell;
    }

    tr.appendChild(td(row.trade_name));
    tr.appendChild(td(row.ingredient));
    tr.appendChild(td(row.route));
    tr.appendChild(td(row.dose_form));
    tr.appendChild(td(row.strength));

    const appParts = [];
    if (row.appl_type || row.appl_no) {
      appParts.push([row.appl_type, row.appl_no].filter(Boolean).join(" "));
    }
    if (row.product_no) {
      appParts.push("Prod " + row.product_no);
    }
    const appStr = appParts.join(" / ");

    tr.appendChild(td(appStr || ""));
    tr.appendChild(td(row.approval_date || ""));

    const marketingParts = [];
    if (row.type) marketingParts.push(row.type);
    if (row.te_code) marketingParts.push("TE " + row.te_code);
    if (row.rld) marketingParts.push("RLD");
    const marketing = marketingParts.join(" / ");

    tr.appendChild(td(marketing || ""));
    prodBody.appendChild(tr);
  });

  return {
    approvalsCount: approvals.length,
    productsCount: obProducts.length
  };
}

// Indications table
function buildIndicationsTable(tbody, data) {
  tbody.innerHTML = "";
  if (!data) return 0;
  const rows = data.indications || [];

  // HTML header:
  //   Disease / Condition | Relation (OMOP) | SNOMED ID | UMLS CUI | Semantic Type
  rows.forEach(function (r) {
    const tr = document.createElement("tr");
    function td(text) {
      const cell = document.createElement("td");
      cell.textContent = text || "-";
      return cell;
    }

    const disease =
      r.concept_name ||
      r.snomed_full_name ||
      "";

    const relation = r.relationship_name || "indication";

    const snomedId =
      r.snomed_conceptid != null && r.snomed_conceptid !== 0
        ? String(r.snomed_conceptid)
        : "";

    const umlsCui = r.umls_cui || "";

    const sty = r.cui_semantic_type || "";

    tr.appendChild(td(disease));
    tr.appendChild(td(relation));
    tr.appendChild(td(snomedId));
    tr.appendChild(td(umlsCui));
    tr.appendChild(td(sty));
    tbody.appendChild(tr);
  });
  return rows.length;
}

// Targets / bioactivity table
function buildTargetsTable(tbody, data) {
  tbody.innerHTML = "";
  if (!data) return 0;
  const rows = data.targets || [];

  // HTML header:
  //   Target (Name / Gene) | Target Class | Pharos (TDL) | UniProt (Accession / SwissProt) |
  //   Action | Activity Type | Activity Value (with Unit) | Mechanism (MoA vs Off-target)
  rows.forEach(function (r) {
    const tr = document.createElement("tr");
    function td(text) {
      const cell = document.createElement("td");
      cell.textContent = text || "-";
      return cell;
    }

    const targetName = r.target_name || r.name || "";
    const gene = r.gene || r.gene_symbol || "";
    let targetDisplay = "";
    if (targetName && gene && gene !== targetName) {
      targetDisplay = targetName + " (" + gene + ")";
    } else if (targetName) {
      targetDisplay = targetName;
    } else {
      targetDisplay = gene || "";
    }

    const tclass = r.target_class || r.class_name || "";
    const tdl = r.tdl || "";
    const acc = r.accession || "";
    const swiss = r.swissprot || "";
    const accDisplay = [acc, swiss].filter(Boolean).join(" / ");

    const action = r.action_type || "";
    const actType = r.act_type || "";

    let actValueStr = "";
    if (r.act_value != null && r.act_value !== "") {
      actValueStr = String(r.act_value);
      if (r.act_unit) actValueStr += " " + r.act_unit;
      if (r.relation) actValueStr = r.relation + " " + actValueStr;
    }

    let mechanism = "";
    if (r.is_moa) {
      mechanism = "MoA (primary mechanism)";
    } else {
      mechanism = "Non-MoA / off-target";
    }

    tr.appendChild(td(targetDisplay));
    tr.appendChild(td(tclass));
    tr.appendChild(td(tdl));
    tr.appendChild(td(accDisplay));
    tr.appendChild(td(action));
    tr.appendChild(td(actType));
    tr.appendChild(td(actValueStr));
    tr.appendChild(td(mechanism));
    tbody.appendChild(tr);
  });
  return rows.length;
}

function buildPharmTable(tbody, data) {
  tbody.innerHTML = "";
  if (!data) return 0;
  const rows = data.pharm_actions || [];
  rows.forEach(function (r) {
    const tr = document.createElement("tr");
    function td(text) {
      const cell = document.createElement("td");
      cell.textContent = text || "-";
      return cell;
    }
    tr.appendChild(td(r.type));
    tr.appendChild(td(r.name));
    tr.appendChild(td(r.class_code));
    tr.appendChild(td(r.source));
    tbody.appendChild(tr);
  });
  return rows.length;
}

function buildSafetySection(container, faersBody, data) {
  container.innerHTML = "";
  faersBody.innerHTML = "";
  if (!data) return { warningsCount: 0, faersCount: 0 };

  const safety = data.safety || {};
  const boxed = safety.boxed_warnings || [];
  const contraindications = safety.contraindications || [];
  const notes = safety.notes || [];

  const faersRows = [];
  (safety.faers_signals || []).forEach(function (r) {
    faersRows.push(Object.assign({}, r));
  });
  (safety.faers_female || []).forEach(function (r) {
    const row = Object.assign({}, r, { sex: "F" });
    faersRows.push(row);
  });
  (safety.faers_male || []).forEach(function (r) {
    const row = Object.assign({}, r, { sex: "M" });
    faersRows.push(row);
  });
  (safety.faers_pediatric || []).forEach(function (r) {
    const row = Object.assign({}, r, { pediatric_flag: true });
    faersRows.push(row);
  });
  (safety.faers_geriatrics || []).forEach(function (r) {
    const row = Object.assign({}, r, { sex: "geriatrics" });
    faersRows.push(row);
  });

  boxed.forEach(function (bw) {
    const div = document.createElement("div");
    div.className = "bb-warning-box";
    const header = document.createElement("div");
    header.className = "bb-warning-header";
    header.textContent = bw.section_title || bw.title || "Black Box Warning";
    const text = document.createElement("div");
    text.className = "bb-warning-text";
    text.textContent = bw.text || "(No summary text available.)";
    div.appendChild(header);
    div.appendChild(text);
    if (bw.pdf_url) {
      const link = document.createElement("div");
      link.className = "bb-warning-link";
      link.innerHTML =
        '<a href="' +
        esc(bw.pdf_url) +
        '" target="_blank" rel="noopener">View FDA Label PDF</a>';
      div.appendChild(link);
    }
    container.appendChild(div);
  });

  contraindications.forEach(function (c) {
    const div = document.createElement("div");
    div.className = "info-panel";
    const header = document.createElement("div");
    header.className = "bb-warning-header";
    header.textContent = c.title || "Contraindication";
    const text = document.createElement("div");
    text.className = "bb-warning-text";
    text.textContent = c.description || c.text || "";
    div.appendChild(header);
    div.appendChild(text);
    container.appendChild(div);
  });

  notes.forEach(function (n) {
    const div = document.createElement("div");
    div.className = "info-panel";
    div.textContent = n;
    container.appendChild(div);
  });

  faersRows.forEach(function (r) {
    const tr = document.createElement("tr");
    function td(text) {
      const cell = document.createElement("td");
      cell.textContent = text || "-";
      return cell;
    }
    tr.appendChild(td(r.level));
    tr.appendChild(td(r.meddra_name));
    tr.appendChild(td(r.meddra_code));
    tr.appendChild(td(r.llr != null ? String(r.llr) : ""));
    tr.appendChild(td(r.llr_threshold != null ? String(r.llr_threshold) : ""));
    const sexAge = [];
    if (r.sex) sexAge.push(r.sex);
    if (r.pediatric_flag) sexAge.push("pediatric");
    tr.appendChild(td(sexAge.join(" / ")));
    faersBody.appendChild(tr);
  });

  return {
    warningsCount: boxed.length + contraindications.length,
    faersCount: faersRows.length
  };
}

function buildCrossrefTable(tbody, data) {
  tbody.innerHTML = "";
  if (!data) return 0;
  const cr = data.cross_refs || {};
  const rows = [];

  if (Array.isArray(cr.identifiers)) {
    cr.identifiers.forEach(function (id) {
      if (!id || !id.identifier) return;
      rows.push({
        label: mapIdType(id.id_type),
        value: String(id.identifier)
      });
    });
  }

  if (cr.atc_codes && cr.atc_codes.length) {
    rows.push({
      label: "ATC Codes",
      value: cr.atc_codes.join(", ")
    });
  }

  if (cr.drug_classes && cr.drug_classes.length) {
    rows.push({
      label: "Drug classes",
      value: cr.drug_classes.join(", ")
    });
  }

  rows.forEach(function (r) {
    const tr = document.createElement("tr");
    const td1 = document.createElement("td");
    const td2 = document.createElement("td");
    td1.textContent = r.label;
    td2.textContent = r.value;
    tr.appendChild(td1);
    tr.appendChild(td2);
    tbody.appendChild(tr);
  });

  return rows.length;
}

// =============================
// Populate a selected drug
// =============================
async function populateDrug(drugKey) {
  if (typeof drug_centralData !== "object" || drug_centralData === null) return;
  if (!drug_centralData[drugKey]) return;

  currentDrug = drugKey;
  const data = drug_centralData[drugKey];

  // Overview
  (function () {
    const section = document.getElementById("overviewSection");
    const container = document.getElementById("overviewCard");
    const headerTitle = section.querySelector(".collapsible-header .title");
    const countEl = section.querySelector(".collapsible-header .count");
    const dg = data.drug || {};
    const hasAnything =
      dg &&
      (dg.name ||
        dg.cas_number ||
        dg.smiles ||
        dg.struct_id != null);

    if (!hasAnything) {
      section.style.display = "none";
      return;
    }
    headerTitle.textContent = "Overview and Structure";
    countEl.textContent = "";
    section.style.display = "block";
    buildOverviewCard(container, data, drugKey);
  })();

  // Regulatory + products
  (function () {
    const section = document.getElementById("regulatorySection");
    const headerTitle = section.querySelector(".collapsible-header .title");
    const countEl = section.querySelector(".collapsible-header .count");
    const regBody = document.getElementById("regulatoryTable");
    const prodBody = document.getElementById("productsTable");
    const reg = data.regulatory || {};
    const hasData =
      (reg.approvals && reg.approvals.length) ||
      (reg.orange_book_products && reg.orange_book_products.length);

    if (!hasData) {
      section.style.display = "none";
      return;
    }
    headerTitle.textContent = "Regulatory Status and Pharmaceutical Products";
    const counts = buildRegulatoryTables(regBody, prodBody, data);
    const totalCount = (counts.approvalsCount || 0) + (counts.productsCount || 0);
    countEl.textContent = String(totalCount);
    section.style.display = "block";
  })();

  // Indications / Drug Use
  (function () {
    const section = document.getElementById("indicationsSection");
    const headerTitle = section.querySelector(".collapsible-header .title");
    const countEl = section.querySelector(".collapsible-header .count");
    const tbody = document.getElementById("indicationsTable");
    const inds = data.indications || [];
    if (!inds.length) {
      section.style.display = "none";
      return;
    }
    headerTitle.textContent = "Drug Use and Indications";
    const n = buildIndicationsTable(tbody, data);
    countEl.textContent = String(n);
    section.style.display = "block";
  })();

  // Targets, bioactivity, pharmacologic actions
  (function () {
    const section = document.getElementById("targetsSection");
    const headerTitle = section.querySelector(".collapsible-header .title");
    const countEl = section.querySelector(".collapsible-header .count");
    const tBody = document.getElementById("targetsTable");
    const pBody = document.getElementById("pharmTable");
    const tRows = data.targets || [];
    const pRows = data.pharm_actions || [];
    if (!tRows.length && !pRows.length) {
      section.style.display = "none";
      return;
    }
    headerTitle.textContent = "Targets and Bioactivity";
    const nt = buildTargetsTable(tBody, data);
    const np = buildPharmTable(pBody, data);
    countEl.textContent = String((nt || 0) + (np || 0));
    section.style.display = "block";
  })();

  // Safety
  (function () {
    const section = document.getElementById("safetySection");
    const headerTitle = section.querySelector(".collapsible-header .title");
    const countEl = section.querySelector(".collapsible-header .count");
    const container = document.getElementById("boxedWarningsContainer");
    const faersBody = document.getElementById("faersTable");
    const safety = data.safety || {};
    const hasSafety =
      (safety.boxed_warnings && safety.boxed_warnings.length) ||
      (safety.contraindications && safety.contraindications.length) ||
      (safety.faers_signals && safety.faers_signals.length) ||
      (safety.faers_female && safety.faers_female.length) ||
      (safety.faers_male && safety.faers_male.length) ||
      (safety.faers_pediatric && safety.faers_pediatric.length) ||
      (safety.faers_geriatrics && safety.faers_geriatrics.length) ||
      (safety.notes && safety.notes.length);

    if (!hasSafety) {
      section.style.display = "none";
      return;
    }
    headerTitle.textContent = "Safety and Adverse Events";
    const counts = buildSafetySection(container, faersBody, data);
    const total = (counts.warningsCount || 0) + (counts.faersCount || 0);
    countEl.textContent = String(total);
    section.style.display = "block";
  })();

  // Cross references
  (function () {
    const section = document.getElementById("crossrefSection");
    const headerTitle = section.querySelector(".collapsible-header .title");
    const countEl = section.querySelector(".collapsible-header .count");
    const tbody = document.getElementById("crossrefTable");
    const cr = data.cross_refs || {};
    const hasCr =
      (cr.identifiers && cr.identifiers.length) ||
      (cr.atc_codes && cr.atc_codes.length) ||
      (cr.drug_classes && cr.drug_classes.length);
    if (!hasCr) {
      section.style.display = "none";
      return;
    }
    headerTitle.textContent = "Cross References and Drug Classes";
    const n = buildCrossrefTable(tbody, data);
    countEl.textContent = String(n);
    section.style.display = "block";
  })();
}

// =============================
// DOM Ready + Theme Sync
// =============================
document.addEventListener("DOMContentLoaded", async function () {
  try {
    RDKit = await initRDKitModule();
  } catch (e) {
    console.warn("RDKit init failed:", e);
  }

  const drugSelect = document.getElementById("drugSelect");
  const validDrugs = [];

  let dataObj = null;
  if (typeof drug_centralData === "object" && drug_centralData !== null) {
    dataObj = drug_centralData;
  }

  if (dataObj) {
    Object.keys(dataObj).forEach(function (key) {
      const d = dataObj[key];
      if (!d || typeof d !== "object") return;
      const dg = d.drug || {};
      const anyEvidence =
        (dg && (dg.name || dg.cas_number || dg.smiles || dg.struct_id != null)) ||
        (d.indications && d.indications.length) ||
        (d.targets && d.targets.length) ||
        (d.regulatory &&
          ((d.regulatory.approvals || []).length ||
           (d.regulatory.orange_book_products || []).length));
      if (anyEvidence) {
        const opt = document.createElement("option");
        opt.value = key;
        const labelName = dg.name || key;
        opt.textContent = labelName;
        drugSelect.appendChild(opt);
        validDrugs.push(key);
      }
    });
  }

  if (validDrugs.length > 0) {
    drugSelect.value = validDrugs[0];
    populateDrug(validDrugs[0]);
  } else {
    drugSelect.style.display = "none";
    document.querySelectorAll(".collapsible-section").forEach(function (sec) {
      sec.style.display = "none";
    });
  }

  // Collapsible behavior
  document.querySelectorAll(".collapsible-header").forEach(function (header) {
    const section = header.closest(".collapsible-section");
    const content = section.querySelector(".collapsible-content");
    const icon = header.querySelector(".toggle-icon");

    header.addEventListener("click", function () {
      const isOpen = content.classList.toggle("open");
      header.classList.toggle("active", isOpen);
      icon.textContent = isOpen ?  "▾" : "▸";
      content.style.maxHeight = "";
    });

    // initial state: closed
    content.classList.remove("open");
    content.style.maxHeight = "";
    icon.textContent = "▸";
  });

  drugSelect.addEventListener("change", function (e) {
    populateDrug(e.target.value);
  });

  // Theme sync via postMessage
  window.addEventListener("message", function (e) {
    if (!e.data || typeof e.data.theme !== "string") return;
    if (e.data.theme === "light") {
      document.body.classList.add("light");
    } else {
      document.body.classList.remove("light");
    }
  });
});
