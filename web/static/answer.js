// static/answer.js
// Author: Ksenia Sokolova
// Updated: 2025-10-21
// static/answer.js


// --- pull statements robustly ---
function getStatements(st){
    const ver = st?.verification || {};
    if (Array.isArray(ver.statements) && ver.statements.length) return ver.statements;
  
    const pairs = ver.pairs || ver.items || ver.claims || [];
    if (Array.isArray(pairs) && pairs.length){
      return pairs.map((it) => {
        if (Array.isArray(it)) {
          const [proof, text] = it;
          const proofs = proof ? (Array.isArray(proof) ? proof : [proof]) : [];
          return { text, proofs };
        }
        const text = it?.text ?? "";
        const pr   = it?.proofs ?? (it?.proof ? [it.proof] : []);
        const proofs = Array.isArray(pr) ? pr : [];
        return { ...it, text, proofs }; // <-- spread first, then canonical fields win
      });
    }
    return [];
  }
  

  
export function renderAnswerVerified(st) {
    const card = byId("answerCard");
    const body = byId("answer");
    if (!card || !body) return;
  
    ensureStyles();
    ensureModal();          // proofs modal
    ensureVerifModal();     // verification summary modal
  
    const ver = st?.verification || {};
    const verdict = String(ver.verdict || "pass").toLowerCase();
    const statements = getStatements(st);

  
    // --- Answer topbar ---
    let top = card.querySelector(".answer-topbar");
    if (!top) {
      top = document.createElement("div");
      top.className = "answer-topbar";
      card.insertBefore(top, body);
    }
    const badgeClass =
      verdict === "pass"    ? "badge badge-supported"  :
      verdict === "partial" ? "badge badge-partial"    :
                               "badge badge-unsupported";
  
    const hasOverall = !!ver?.evidence?.overall || !!ver?.overall;
    top.innerHTML = `
      <div class="row" style="gap:8px; align-items:center; flex-wrap:wrap;">
        <span class="${badgeClass}">VERDICT: ${verdict.toUpperCase()}</span>
        <span class="muted">Hover for feedback; click a highlight to see references to data.</span>
        ${hasOverall ? `<button class="ghost small" id="verifDetailsBtn" aria-haspopup="dialog">Verification details</button>` : ""}
        <label class="toggle small" style="margin-left:auto;">
          <input type="checkbox" id="hlToggle" />
          <span id="hlToggleLabel">Show highlights</span>
        </label>
      </div>
    `;
  
 // --- build body (statement-based, no manual <br>) ---
const HIGHLIGHT_MIN_WORDS = 3;
body.innerHTML = "";
const frag = document.createDocumentFragment();

const isBlockyMdStart = (txt) => /^(#{1,6}\s+|[-*]\s|\d+\.\s|>|\s*```|={3,}\s*$|-{3,}\s*$)/m.test(txt.trimStart());

// track previous inline chip so we can glue punctuation & preserve boundary spaces
let prevSpan = null;
let prevLiteral = "";

statements.forEach((s, idx) => {
    const raw = normalizeRawText(s?.text || "");
    let literal = raw; // DO NOT trim/collapse
  
    // full-block chunks (start with MD block mark)
    const isBlockyMdStart = (t) => /^(#{1,6}\s+|[-*]\s|\d+\.\s|>|\s*```)/.test(t);
  
  
    // A) Whole chunk is a block → render as block & reset inline context.
    if (isBlockyMdStart(literal)) {
      const block = document.createElement("div");
      block.className = "answer-block-md";
      block.innerHTML = mdToHTML(literal);
      frag.appendChild(block);
      prevSpan = null;
      prevLiteral = "";
      return;
    }
  
    // B) Chunk contains a block start later, e.g. ".\n\n## Heading..."
    const blockBoundary = literal.search(/\n{2,}(?=(#{1,6}\s+|[-*]\s|\d+\.\s|>|\s*```))/);
    if (blockBoundary >= 0) {
      const prefix = literal.slice(0, blockBoundary);      // inline before the block (maybe just ".")
      const rest   = literal.slice(blockBoundary + 2);     // starts with the block
  
      // 1) If prefix is punctuation-only, merge into previous chip to avoid wraps/gaps.
      if (/^[\s]*[.,:;!?)]\s*$/.test(prefix) && prevSpan) {
        const patched = (/^ /.test(prefix) && !/\s$/.test(prevLiteral))
          ? "\u00A0" + prefix.slice(1)   // preserve a single intended leading space
          : prefix;
        prevSpan.innerHTML = prevSpan.innerHTML + mdInlineHTML(patched);
        prevLiteral += patched;
      } else if (prefix) {
        // render non-empty non-punct prefix as a normal inline chip
        const span = document.createElement("span");
        span.className = "answer-chunk";
        span.dataset.idx = String(idx);
        span.dataset.statement = JSON.stringify(s);
  
        const llm = String(s?.llm_label || "").toLowerCase();
        const v   = String(s?.verdict   || "").toLowerCase();
        let colorClass = "stmt-unsupported";
        if (llm === "unsupported") colorClass = "stmt-llm-unsupported";
        else if (v === "supported") colorClass = "stmt-supported";
        else if (v === "partial")   colorClass = "stmt-partial";
        if (wordCount(prefix) >= HIGHLIGHT_MIN_WORDS) span.classList.add("hl-chip", colorClass);
        else span.classList.add("hl-ghost");
  
        const prefixRendered = (/^ /.test(prefix) && prevSpan && !/\s$/.test(prevLiteral))
          ? "\u00A0" + prefix.slice(1)
          : prefix;
        span.innerHTML = mdInlineHTML(prefixRendered);
        frag.appendChild(span);
        prevSpan = span;
        prevLiteral = prefixRendered;
      }
  
      // 2) From `rest`, extract just the FIRST block element (e.g., the heading line)
      //    and keep the remaining tail as inline (so it can continue with next chunks).
      const mHeading = rest.match(/^(#{1,6}\s+[^\n]+)(\n+)?/);
      if (mHeading) {
        const headingMd   = mHeading[1];           // "## Mechanistic Rationale"
        let tailInline    = rest.slice(mHeading[0].length); // whatever follows heading line breaks
  
        // Render the heading as a block
        const block = document.createElement("div");
        block.className = "answer-block-md";
        block.innerHTML = mdToHTML(headingMd);
        frag.appendChild(block);
  
        // Strip ONLY the leading newlines that belonged to MD separation after the heading.
        tailInline = tailInline.replace(/^\n+/, "");
  
        // If there is tail text (e.g., "TP53's involvement in"), render it as inline
        if (tailInline) {
          const span = document.createElement("span");
          span.className = "answer-chunk";
          span.dataset.idx = String(idx);
          span.dataset.statement = JSON.stringify(s);
  
          const llm = String(s?.llm_label || "").toLowerCase();
          const v   = String(s?.verdict   || "").toLowerCase();
          let colorClass = "stmt-unsupported";
          if (llm === "unsupported") colorClass = "stmt-llm-unsupported";
          else if (v === "supported") colorClass = "stmt-supported";
          else if (v === "partial")   colorClass = "stmt-partial";
          if (wordCount(tailInline) >= HIGHLIGHT_MIN_WORDS) span.classList.add("hl-chip", colorClass);
          else span.classList.add("hl-ghost");
  
          // If the tail begins with one real space and the previous node is a block,
          // HTML would otherwise collapse it at the boundary -> NBSP just for the first one.
          if (/^ /.test(tailInline) && !prevSpan) {
            tailInline = "\u00A0" + tailInline.slice(1);
          }
  
          span.innerHTML = mdInlineHTML(tailInline);
          frag.appendChild(span);
          prevSpan = span;
          prevLiteral = tailInline;
        } else {
          // heading ended the chunk; reset inline context
          prevSpan = null;
          prevLiteral = "";
        }
        return;
      }
  
      // If it's another kind of block (list/code/blockquote), fall back to full md block:
      const block = document.createElement("div");
      block.className = "answer-block-md";
      block.innerHTML = mdToHTML(rest);
      frag.appendChild(block);
      prevSpan = null;
      prevLiteral = "";
      return;
    }
  
    // C) No block inside — pure inline chip handling
    if (/^[\s]*[.,:;!?)]\s*$/.test(literal) && prevSpan) {
      const patched = (/^ /.test(literal) && !/\s$/.test(prevLiteral))
        ? "\u00A0" + literal.slice(1)
        : literal;
      prevSpan.innerHTML = prevSpan.innerHTML + mdInlineHTML(patched);
      prevLiteral += patched;
      return;
    }
  
    if (/^ /.test(literal) && prevSpan && !/\s$/.test(prevLiteral)) {
      literal = "\u00A0" + literal.slice(1);
    }
  
    const span = document.createElement("span");
    span.className = "answer-chunk";
    span.dataset.idx = String(idx);
    span.dataset.statement = JSON.stringify(s);
  
    const llm = String(s?.llm_label || "").toLowerCase();
    const v   = String(s?.verdict   || "").toLowerCase();
    let colorClass = "stmt-unsupported";
    if (llm === "unsupported") colorClass = "stmt-llm-unsupported";
    else if (v === "supported") colorClass = "stmt-supported";
    else if (v === "partial")   colorClass = "stmt-partial";
    if (wordCount(literal) >= HIGHLIGHT_MIN_WORDS) span.classList.add("hl-chip", colorClass);
    else span.classList.add("hl-ghost");
  
    span.innerHTML = mdInlineHTML(literal);
    frag.appendChild(span);
    prevSpan = span;
    prevLiteral = literal;
  });
  

body.appendChild(frag);
card.style.display = statements.length ? "block" : "none";




    // interactions
    wireTooltips(body);
  
    body.querySelectorAll(".answer-chunk").forEach((el) => {
      el.addEventListener("click", () => {
        if (card.dataset.highlights !== "on") return; // click disabled when OFF
        const s = safeJson(el.dataset.statement) || {};
        openProofsModal(s);
      });
    });
  
    // Verification details
    const btn = byId("verifDetailsBtn");
    if (btn) btn.addEventListener("click", () => openVerifModal(ver));
  
    // Highlight toggle (default ON)
const tgl = byId("hlToggle");
const lbl = byId("hlToggleLabel");
card.dataset.highlights = "on";
if (tgl && lbl) {
  tgl.checked = true; // ON by default

  const updateLabel = () => {
    lbl.textContent = tgl.checked ? "Hide verifier highlights" : "Show verifier highlights";
  };
  updateLabel();

  tgl.addEventListener("change", () => {
    card.dataset.highlights = tgl.checked ? "on" : "off";
    updateLabel();
    if (!tgl.checked && typeof TIP !== "undefined" && TIP) {
      TIP.classList.remove("show");
      TIP.style.left = TIP.style.top = "";
    }
  });
}
  }
  
  /* ---------------- utils ---------------- */
  function byId(id){ return document.getElementById(id); }
  function esc(s){
    return String(s ?? "").replace(/[&<>"']/g, ch => (
      {'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[ch]
    ));
  }
  function wordCount(s){ return s.trim().split(/\s+/).filter(Boolean).length; }
  function safeJson(s){ try{ return JSON.parse(String(s||"{}")); }catch{ return null; } }
  function makeTooltipHTML({ feedback, reasons, speculation }){
    const pills = (reasons||[]).map(r => `<span class="pill">${esc(String(r))}</span>`).join(" ");
    const spec  = speculation ? `<span class="pill">speculation</span>` : "";
    return `
      <div class="tip-line">${esc(feedback || "No issues detected")}</div>
      <div class="tip-line pills">${spec}${pills ? " " + pills : ""}</div>
    `;
  }
  
  /* =====================================================================
     MARKDOWN: replaced custom regex with markdown-it + DOMPurify
     ===================================================================== */
  
  /* Singleton markdown-it instance configured for our classes + safe links */
  function getMd() {
    if (window.__MD__) return window.__MD__;
  
    const md = window.markdownit({
      html: false,       // disallow raw HTML from model text
      linkify: true,     // autolink URLs
      typographer: true, // smart quotes/dashes
      breaks: false
    });
  
    // Ensure all links open safely in a new tab
    const defOpen = md.renderer.rules.link_open || ((t, i, o, e, s) => s.renderToken(t, i, o));
    md.renderer.rules.link_open = function(tokens, idx, options, env, self) {
      const a = tokens[idx];
      const tIdx = a.attrIndex('target');
      const rIdx = a.attrIndex('rel');
      if (tIdx < 0) a.attrPush(['target','_blank']); else a.attrs[tIdx][1] = '_blank';
      if (rIdx < 0) a.attrPush(['rel','noopener noreferrer']); else a.attrs[rIdx][1] = 'noopener noreferrer';
      return defOpen(tokens, idx, options, env, self);
    };
  
    // Add our heading classes
    md.renderer.rules.heading_open = function(tokens, idx, options, env, self) {
      const token = tokens[idx];
      const level = token.tag.replace('h',''); // "1".."6"
      token.attrJoin('class', `md-h${level}`);
      return self.renderToken(tokens, idx, options);
    };
  
    // Inline/Block code classes
    md.renderer.rules.code_inline = function(tokens, idx) {
      return `<code class="md-code-inline">${escapeHtml(tokens[idx].content)}</code>`;
    };
    md.renderer.rules.code_block = function(tokens, idx) {
      return `<pre class="md-code"><code>${escapeHtml(tokens[idx].content)}</code></pre>`;
    };
    md.renderer.rules.fence = function(tokens, idx) {
      const token = tokens[idx];
      const info = (token.info || '').trim();
      const langClass = info ? ` class="language-${escapeHtml(info)}"` : '';
      return `<pre class="md-code"><code${langClass}>${escapeHtml(token.content)}</code></pre>`;
    };
  
    window.__MD__ = md;
    return md;
  }
  
  function escapeHtml(s) {
    return String(s ?? "").replace(/[&<>"']/g, ch => (
      {'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[ch]
    ));
  }
  
  function sanitize(html, { inline = false } = {}) {
    const ALLOWED_TAGS = inline
      ? ['em','strong','code','a','del','span','sub','sup','br']
      : ['p','em','strong','code','pre','a','ul','ol','li','h1','h2','h3','h4','h5','h6','blockquote','del','span','sub','sup','br'];
    return DOMPurify.sanitize(html, {
      ALLOWED_TAGS,
      ALLOWED_ATTR: ['href','title','target','rel','class']
    });
  }
  
  /* Inline-only Markdown for chips */
  function mdInlineHTML(src){
    const html = getMd().renderInline(String(src ?? ''));
    return sanitize(html, { inline: true }); 
  }
  
  /* Block Markdown outside chips */
  function mdToHTML(src){
    const html = getMd().render(String(src ?? ''));
    return sanitize(html).trim();
  }
  
function normalizeRawText(src) {
    return String(src ?? "").replace(/\r\n?/g, "\n");
  }
  
  
  /* ---------------- tooltip ---------------- */
  let TIP;
  function wireTooltips(scope){
    if (!TIP){
      TIP = document.createElement("div");
      TIP.className = "cite-tip";
      document.body.appendChild(TIP);
    }
    scope.querySelectorAll(".answer-chunk").forEach(node => {
      node.addEventListener("mouseenter", (e) => {
        if (byId("answerCard")?.dataset.highlights !== "on") return;
        // Generate tooltip HTML dynamically from statement data to avoid HTML escaping issues
        const stmt = safeJson(node.dataset.statement || "{}");
        if (!stmt) return;
        const html = makeTooltipHTML({
          feedback: stmt.llm_feedback,
          reasons: Array.isArray(stmt.reasons) ? stmt.reasons : [],
          speculation: !!stmt.is_speculation
        });
        if (!html) return;
        TIP.innerHTML = html;
        TIP.classList.add("show");
        positionTip(e);
      });
      node.addEventListener("mousemove", (e) => {
        if (byId("answerCard")?.dataset.highlights !== "on") return;
        positionTip(e);
      });
      node.addEventListener("mouseleave", () => {
        TIP.classList.remove("show"); TIP.style.left = TIP.style.top = "";
      });
    });
  }
  function positionTip(e){
    const pad = 12;
    TIP.style.left = (e.clientX + pad) + "px";
    TIP.style.top  = (e.clientY + pad) + "px";
  }
  
  /* ---------------- proofs modal ---------------- */
  let MODAL;
  function ensureModal(){
    if (MODAL) return;
    MODAL = document.createElement("div");
    MODAL.id = "proofsModal";
    MODAL.innerHTML = `
      <div class="modal-backdrop"></div>
      <div class="modal-sheet">
        <div class="modal-head">
          <div class="modal-title">References</div>
          <button class="ghost modal-close" aria-label="Close">✕</button>
        </div>
        <div class="modal-body"></div>
      </div>
    `;
    document.body.appendChild(MODAL);
    MODAL.querySelector(".modal-backdrop").addEventListener("click", closeModal);
    MODAL.querySelector(".modal-close").addEventListener("click", closeModal);
    document.addEventListener("keydown", (e) => e.key === "Escape" && closeModal());
  }
  function openProofsModal(stmt){
    if (!stmt) return;
    ensureModal();  // Initialize modal if not already created
    const body = MODAL.querySelector(".modal-body");
    const proofs = Array.isArray(stmt?.proofs) ? stmt.proofs : [];
    const reasons = Array.isArray(stmt?.reasons) ? stmt.reasons : [];
    const feedback = stmt?.llm_feedback;
  
    const list = proofs.map((p, i) => {
      const title = (p?.title || p?.source || `Citation ${i+1}`);
      const cited = p?.cited_text
        ? `<blockquote class="cite-block">${esc(String(p.cited_text))}</blockquote>`
        : `<div class="muted">No cited text provided.</div>`;
      const url = p?.url ? `<div class="muted sm"><a href="${esc(p.url)}" target="_blank" rel="noopener">Open source ↗</a></div>` : "";
      return `
        <li class="ref-item">
          <details class="ref-details" ${i===0 ? "open" : ""}>
            <summary class="ref-summary">${esc(title)}</summary>
            ${cited}
            ${url}
          </details>
        </li>
      `;
    }).join("");
  
    const reasonsPills = reasons.map(r => `<span class="pill">${esc(String(r))}</span>`).join(" ");
    const spec = stmt?.is_speculation ? `<span class="pill">speculation</span>` : "";
  
    body.innerHTML = `
      <div class="muted" style="margin-bottom:6px;">Statement</div>
      <div class="stmt-text">${mdToHTML(stmt?.text || "")}</div>
  
      <div class="muted" style="margin:10px 0 6px;">Citations</div>
      <ol class="refs-list">${list || `<li class="muted">No proofs attached.</li>`}</ol>
  
      <div class="muted" style="margin:10px 0 6px;">Feedback & reasons</div>
      <div class="mono wrap">${feedback ? mdToHTML(feedback) : "—"}</div>
      <div>${spec} ${reasonsPills}</div>
    `;
    MODAL.classList.add("open");
  }
  function closeModal(){ MODAL?.classList.remove("open"); }
  
  /* ---------------- verification summary modal ---------------- */
  let VERIF_MODAL;
  function ensureVerifModal(){
    if (VERIF_MODAL) return;
    VERIF_MODAL = document.createElement("div");
    VERIF_MODAL.id = "verifModal";
    VERIF_MODAL.innerHTML = `
      <div class="modal-backdrop"></div>
      <div class="modal-sheet">
        <div class="modal-head">
          <div class="modal-title">Verification details</div>
          <button class="ghost modal-close" aria-label="Close">✕</button>
        </div>
        <div class="modal-body"></div>
      </div>
    `;
    document.body.appendChild(VERIF_MODAL);
    VERIF_MODAL.querySelector(".modal-backdrop").addEventListener("click", () => closeVerif());
    VERIF_MODAL.querySelector(".modal-close").addEventListener("click", () => closeVerif());
    document.addEventListener("keydown", (e) => e.key === "Escape" && closeVerif());
  }
  function openVerifModal(ver){
    if (!ver) return;
    ensureVerifModal();  // Initialize modal if not already created
    const body = VERIF_MODAL.querySelector(".modal-body");
    const overall = ver?.evidence?.overall || ver?.overall || {};
    const quality = overall?.support_quality ? `<span class="pill">${esc(overall.support_quality)}</span>` : "";
    const summary = overall?.summary ? `<p class="wrap">${esc(overall.summary)}</p>` : `<p class="muted">No overall summary.</p>`;
    const concerns = Array.isArray(overall?.concerns) && overall.concerns.length
      ? `<ul class="md-list">${overall.concerns.map(c => `<li>${esc(String(c))}</li>`).join("")}</ul>`
      : `<div class="muted">No concerns listed.</div>`;
  
    const v = String(ver?.verdict || "").toLowerCase();
    const verdictClass =
      v === "pass"    ? "badge badge-supported"  :
      v === "partial" ? "badge badge-partial"    :
                        "badge badge-unsupported";
  
    body.innerHTML = `
      <section class="verif-section">
        <h3 class="md-h3">Overall quality: ${quality}</h3>
        ${summary}
      </section>
      <section class="verif-section">
        <h4 class="md-h4">Concerns</h4>
        ${concerns}
      </section>
      ${ver?.verdict ? `
      <section class="verif-section">
        <h4 class="md-h4">Verdict</h4>
        <p><span class="${verdictClass}">VERDICT: ${esc(String(ver.verdict).toUpperCase())}</span></p>
      </section>` : ""}
    `;
    VERIF_MODAL.classList.add("open");
  }
  function closeVerif(){ VERIF_MODAL?.classList.remove("open"); }
  
  /* ---------------- styles ---------------- */
  let STYLED=false;
  function ensureStyles(){
    if (STYLED) return; STYLED = true;
    const css = `
  /* Chips */
  .answer .hl-chip {
    display: inline;
    padding: .05em .15em;
    border-radius: .45em;
    box-decoration-break: clone;
    -webkit-box-decoration-break: clone;
  }
  .answer .hl-ghost { padding: 0 .1em; border-radius: .25em; }
  
  /* Colors (brighter dark; cleaner light) */
  .answer .stmt-supported   { background: rgba(16,185,129,.36); }
  .answer .stmt-partial     { background: rgba(148,163,184,.34); }
  .answer .stmt-unsupported { background: rgba(239,68,68,.36); }
  .answer .stmt-llm-unsupported { background: rgba(190,18,60,.42); }
  
  body.light .answer .stmt-supported   { background: rgba(5,150,105,.22); outline: 1px solid rgba(5,150,105,.28); }
  body.light .answer .stmt-partial     { background: rgba(100,116,139,.20); outline: 1px solid rgba(100,116,139,.26); }
  body.light .answer .stmt-unsupported { background: rgba(220,38,38,.22);  outline: 1px solid rgba(220,38,38,.28); }
  body.light .answer .stmt-llm-unsupported { background: rgba(159,18,57,.26); outline: 1px solid rgba(159,18,57,.30); }
  
  /* Markdown blocks rendered inline */
  .answer-block-md {
    margin: 0.65em 0;
    padding: 0.4em 0;
  }
  .answer-block-md h1,
  .answer-block-md h2,
  .answer-block-md h3,
  .answer-block-md h4,
  .answer-block-md h5,
  .answer-block-md h6 {
    margin: 0.4em 0 0.35em;
    font-weight: 700;
  }
  .answer-block-md h3,
  .answer-block-md h4 {
    border-bottom: 1px solid rgba(148,163,184,.25);
    padding-bottom: 4px;
  }
  .answer-block-md ul,
  .answer-block-md ol {
    margin: 0.35em 0 0.35em 1.3em;
  }
  .answer-block-md pre {
    background: rgba(15,23,42,0.85);
    border-radius: 10px;
    padding: 10px 12px;
    overflow: auto;
    font-size: 13px;
    line-height: 1.4;
    border: 1px solid rgba(255,255,255,.08);
  }
  .answer-block-md code {
    background: rgba(148,163,184,.18);
    border-radius: 6px;
    padding: 2px 6px;
    font-size: 13px;
  }
  .answer-block-md blockquote {
    border-left: 3px solid rgba(148,163,184,.35);
    margin: 8px 0;
    padding: 4px 10px;
    color: var(--muted,#94a3b8);
    font-style: italic;
    background: rgba(148,163,184,.12);
    border-radius: 6px;
  }

  /* Inline code & links inside chips */
  .answer .hl-chip a { text-decoration: underline; }
  .answer .hl-chip .md-code-inline {
    background: rgba(148,163,184,.22);
    padding: .05em .3em; border-radius: .35em;
  }
  .answer { line-height: 1.55; }
  .answer p { margin: 0.35em 0; }
  .answer .hl-chip { margin: 0 1px 1px 0; }
  .answer ul, .answer ol { margin-top: 0.35em; margin-bottom: 0.35em; }
  
  /* Pills / badges */
  .answer-topbar { margin: -2px 0 10px; }
  .badge { display:inline-block; padding:2px 8px; border-radius:999px; font-size:11px; border:1px solid currentColor; }
  .badge-supported  { color:#d1fae5; background:rgba(16,185,129,.20); border-color:rgba(5,150,105,.55); }
  .badge-partial    { color:#8a6100; background:rgba(234,179,8,.20);  border-color:rgba(161,98,7,.40); }
  .badge-unsupported{ color:#f8fafc; background:rgba(255,255,255,.16); border-color:rgba(255,255,255,.32); }
  body.light .badge-supported  { color:#065f46; background:rgba(5,150,105,.14); border-color:rgba(5,150,105,.34); }
  body.light .badge-partial    { color:#854d0e; background:rgba(234,179,8,.16);  border-color:rgba(202,138,4,.34); }
  body.light .badge-unsupported{ color:#1f2937; background:rgba(148,163,184,.20); border-color:rgba(51,65,85,.30); }
  
  button.ghost.small { font-size:12px; padding:4px 8px; border-radius:8px; border:1px solid var(--line,#2a3a72); background:transparent; color:inherit; }
  .toggle.small { display:flex; align-items:center; gap:6px; font-size:12px; }
  .toggle.small input { transform: translateY(1px); }
  
  /* Hide highlight visuals when OFF */
  #answerCard[data-highlights="off"] .answer .hl-chip { background: transparent !important; outline: none !important; }
  #answerCard[data-highlights="off"] .answer .answer-chunk { cursor: text; }
  
  /* Tooltip (wraps) */
  .cite-tip {
    position:fixed; z-index:10000; max-width: 560px;
    background:rgba(12, 11, 10, 0.95); color:#f3f0ea;
    border:1px solid rgba(255,180,120,0.25); border-radius:10px;
    padding:8px 10px; font-size:12px; line-height:1.45;
    box-shadow:0 8px 24px rgba(0,0,0,.4); pointer-events:none;
    opacity:0; transform:translateY(4px);
    transition:opacity .08s ease, transform .08s ease;
    white-space: normal; overflow-wrap: anywhere; word-break: break-word;
  }
  .cite-tip.show { opacity:1; transform:translateY(0); }
  .cite-tip .tip-line + .tip-line { margin-top:4px; }
  .cite-tip .pill { margin:2px 4px 0 0; }

  body.light .cite-tip {
    background: rgba(255, 255, 255, 0.98);
    color: #0b0f1a;
    border-color: rgba(12, 18, 36, 0.2);
    box-shadow: 0 8px 24px rgba(0,0,0,.15);
  }
  
  /* Modals */
  #proofsModal, #verifModal { position: fixed; inset: 0; display: none; z-index: 10001; }
  #proofsModal.open, #verifModal.open { display: flex; align-items: center; justify-content: center; padding: 4vh 14px; }
  #proofsModal .modal-backdrop, #verifModal .modal-backdrop { position: fixed; inset: 0; background: rgba(0,0,0,.65); backdrop-filter: blur(2px); }
  #proofsModal .modal-sheet, #verifModal .modal-sheet {
    position: relative; width: min(100%, 1040px); max-height: 84vh;
    background: #0c0b0a; color: var(--fg, #f3f0ea);
    border: 1px solid rgba(255,180,120,0.20); border-radius: 14px;
    box-shadow: 0 18px 46px rgba(0,0,0,.5);
    display: flex; flex-direction: column; overflow: hidden;
  }
  #proofsModal .modal-head, #verifModal .modal-head { display:flex; align-items:center; justify-content:space-between; gap:8px; padding:12px 14px; border-bottom:1px solid rgba(255,180,120,0.12); }
  #proofsModal .modal-title, #verifModal .modal-title { font-weight:700; }
  #proofsModal .modal-body, #verifModal .modal-body { padding:12px 14px; overflow:auto; }
  body.light #proofsModal .modal-sheet, body.light #verifModal .modal-sheet {
    background:#ffffff; color:#0b1020; border-color:#d8dce6; box-shadow: 0 18px 46px rgba(0,0,0,.2);
  }
  
  /* Proofs list */
  .refs-list { margin:0; padding-left: 20px; }
  .ref-item { margin: 4px 0; }
  .ref-details { border:1px solid rgba(255,180,120,0.12); border-radius:8px; padding:6px 8px; background: rgba(255,200,150,0.03); }
  .ref-details[open] { background: rgba(255,200,150,0.06); }
  .ref-summary { cursor:pointer; font-weight:600; list-style:none; }
  .ref-summary::-webkit-details-marker { display:none; }
  .cite-block {
    margin:8px 0 4px 0; padding:8px 10px;
    border-left:3px solid rgba(255,180,120,0.25);
    background: rgba(255,255,255,.03);
    white-space: pre-wrap; overflow-wrap: anywhere; word-break: break-word;
  }
  
  /* Markdown styles */
  .md-h1,.md-h2,.md-h3,.md-h4,.md-h5,.md-h6 { margin: .35em 0 .2em; line-height:1.2; }
  .md-h1{ font-size:1.6rem; } .md-h2{ font-size:1.35rem; }
  .md-h3{ font-size:1.2rem; }  .md-h4{ font-size:1.05rem; }
  .md-code { background: rgba(148,163,184,.18); padding:10px; border-radius:8px; overflow:auto; }
  .md-code-inline { background: rgba(148,163,184,.18); padding: .1em .35em; border-radius:.35em; }
  .md-list { margin: .3em 0 .3em 1.2em; }
  .wrap { white-space: normal; overflow-wrap: anywhere; word-break: break-word; }
  .stmt-text { white-space: normal; overflow-wrap:anywhere; word-break:break-word; }

  .answer .answer-chunk { white-space: pre-wrap; }
  /* Chips: no horizontal margins that fake a "space" */
.answer .hl-chip,
.answer .hl-ghost {
  margin: 0;                /* was: 0 1px 1px 0 */
}

/* Preserve author spacing/newlines in chips */
.answer .answer-chunk { white-space: pre-wrap; }

/* If a chip starts with punctuation, don't add left padding (no fake gap) */
.answer .answer-chunk.punct-start {
  padding-left: 0 !important;
}

/* If the NEXT chip starts with punctuation, trim current chip's right padding */
.answer .answer-chunk.tight-right {
  padding-right: 0 !important;
}

  
  /* Sections in verification modal */
  .verif-section + .verif-section { margin-top: 12px; }
  `;
    const style = document.createElement("style");
    style.textContent = css;
    document.head.appendChild(style);
  }
  
