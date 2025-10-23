// Author - Ksenia Sokolova
// Date - 2025-10-22
// static/waiting_module.js
// Lightweight waiting overlay + polling helper.
// Works without any CSS from your app (injects its own <style>).

export class Wait {
    static attach(opts = {}) { return new Wait(opts); }
  
    constructor({ parent = document.body } = {}) {
      this._buildDOM(parent);
      this._startTs = 0;
      this._tickTimer = null;
      this._onCancel = null;
    }
  
    _buildDOM(parent) {
      // Styles (scoped enough to avoid collisions)
      const css = `
      .alv-wait-backdrop {
        position: fixed; inset: 0;
        background: rgba(5,8,20,0.55);
        backdrop-filter: blur(6px) saturate(130%);
        display: none; z-index: 9999;
      }
      .alv-wait-panel {
        position: absolute; inset: 0; display: grid; place-items: center;
        padding: 16px;
      }
      .alv-wait-card {
        width: min(520px, 90vw);
        border-radius: 16px;
        background: var(--card, #121832);
        color: var(--fg, #e6e8ef);
        box-shadow: 0 18px 40px rgba(0,0,0,.35);
        border: 1px solid var(--line, #1e2a52);
        padding: 20px;
        font-family: var(--sans, system-ui, -apple-system, Segoe UI, Roboto, Helvetica, Arial);
      }
      .alv-wait-hdr { display: flex; align-items: center; gap: 10px; margin-bottom: 10px; }
      .alv-wait-title { font-weight: 700; font-size: 16px; }
      .alv-wait-sub { color: var(--muted, #aab4d6); font-size: 13px; }
      .alv-wait-row { display: flex; align-items: center; gap: 10px; margin-top: 12px; }
      .alv-wait-spinner {
        width: 22px; height: 22px; border-radius: 50%;
        border: 3px solid rgba(255,255,255,0.25);
        border-top-color: rgba(255,255,255,0.95);
        animation: alvspin 0.9s linear infinite;
      }
      @keyframes alvspin { to { transform: rotate(360deg); } }
  
      .alv-wait-bar-wrap {
        margin-top: 12px; width: 100%;
        height: 8px; border-radius: 999px;
        background: rgba(255,255,255,0.08);
        overflow: hidden;
      }
      .alv-wait-bar {
        width: 12%;
        height: 100%;
        background: linear-gradient(90deg, #6ea8ff, #b1c8ff);
        border-radius: 999px;
        animation: alvbar 1.6s ease-in-out infinite;
      }
      @keyframes alvbar {
        0% { transform: translateX(-60%); width: 20%; }
        50% { transform: translateX(40%); width: 60%; }
        100% { transform: translateX(120%); width: 20%; }
      }
  
      .alv-wait-footer { display: flex; justify-content: space-between; align-items: center; margin-top: 14px; }
      .alv-wait-time { font-size: 12px; color: var(--muted, #aab4d6); }
      .alv-wait-actions { display: flex; gap: 8px; }
      .alv-btn {
        padding: 8px 14px; border: none; border-radius: 999px;
        cursor: pointer; font-weight: 600; font-size: 14px;
        background: var(--accent, #2a3a72); color: var(--fg, #e6e8ef);
        box-shadow: 0 4px 12px rgba(0,0,0,.2);
      }
      .alv-btn.ghost { background: transparent; box-shadow: none; border: 1px solid var(--line, #1e2a52); }
      .alv-wait-err { color: #ffb3b3; margin-top: 10px; font-size: 13px; display: none; }
      `;
      const style = document.createElement('style');
      style.textContent = css;
      document.head.appendChild(style);
  
      // DOM
      this.backdrop = document.createElement('div');
      this.backdrop.className = 'alv-wait-backdrop';
      this.backdrop.setAttribute('role', 'dialog');
      this.backdrop.setAttribute('aria-modal', 'true');
      this.backdrop.innerHTML = `
        <div class="alv-wait-panel">
          <div class="alv-wait-card">
            <div class="alv-wait-hdr">
              <div class="alv-wait-spinner" aria-hidden="true"></div>
              <div>
                <div class="alv-wait-title" id="alv_wait_title">Working…</div>
                <div class="alv-wait-sub" id="alv_wait_msg">Submitting your request</div>
              </div>
            </div>
            <div class="alv-wait-bar-wrap" aria-hidden="true"><div class="alv-wait-bar"></div></div>
            <div class="alv-wait-err" id="alv_wait_err"></div>
            <div class="alv-wait-footer">
              <div class="alv-wait-time"><span id="alv_wait_elapsed">0:00</span> elapsed</div>
              <div class="alv-wait-actions">
                <button class="alv-btn ghost" id="alv_wait_cancel">Cancel</button>
              </div>
            </div>
          </div>
        </div>
      `;
      parent.appendChild(this.backdrop);
  
      this.$title   = this.backdrop.querySelector('#alv_wait_title');
      this.$msg     = this.backdrop.querySelector('#alv_wait_msg');
      this.$err     = this.backdrop.querySelector('#alv_wait_err');
      this.$elapsed = this.backdrop.querySelector('#alv_wait_elapsed');
      this.$cancel  = this.backdrop.querySelector('#alv_wait_cancel');
  
      this.$cancel.addEventListener('click', () => this._handleCancel());
      window.addEventListener('keydown', (e) => {
        if (this.isShown() && e.key === 'Escape') this._handleCancel();
      });
    }
  
    isShown() { return this.backdrop.style.display !== 'none'; }
  
    start(title = 'Working…', msg = 'Please wait') {
      this.$err.style.display = 'none';
      this.$err.textContent = '';
      this.setTitle(title);
      this.setMessage(msg);
      this._startTs = Date.now();
      this._tick();
      this._tickTimer = setInterval(() => this._tick(), 1000);
      this.backdrop.style.display = 'block';
    }
  
    setTitle(t) { this.$title.textContent = t || 'Working…'; }
    setMessage(m) { this.$msg.textContent = m || ''; }
    setError(err) {
      if (err) {
        this.$err.textContent = String(err);
        this.$err.style.display = 'block';
      } else {
        this.$err.style.display = 'none';
        this.$err.textContent = '';
      }
    }
  
    onCancel(fn) { this._onCancel = fn; }
    _handleCancel() { if (typeof this._onCancel === 'function') this._onCancel(); this.stop(); }
  
    stop() {
      clearInterval(this._tickTimer); this._tickTimer = null;
      this.backdrop.style.display = 'none';
    }
  
    _tick() {
      const s = Math.max(0, Math.floor((Date.now() - this._startTs) / 1000));
      const mm = String(Math.floor(s / 60)).padStart(1, '0');
      const ss = String(s % 60).padStart(2, '0');
      this.$elapsed.textContent = `${mm}:${ss}`;
    }
  }
