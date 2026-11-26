# 5-Minute Zero-Cost Cloud Deployment Guide (Google Sheets + Netlify / GitHub Pages)

This lets anyone open ONE HTTPS link on their iPhone or Android phone (zero App Store approval, zero account creation), log check-ins and N-of-1 trials offline or online, and stream every anonymized row straight into your Google Sheet or Supabase database.

---

## Option A (Recommended, $0): Google Sheets Collector
1. Open Google Sheets and create a blank spreadsheet named **`TRPM8 Study Data`**.
2. Go to **Extensions → Apps Script**.
3. Delete the placeholder `myFunction` code and paste the full contents of [`google_apps_script.gs`](./google_apps_script.gs).
4. Click **Deploy → New deployment → Select type: Web app**:
   - **Execute as:** `Me`
   - **Who has access:** `Anyone`
5. Authorize when prompted and copy the **Web app URL** (`https://script.google.com/macros/s/.../exec`).
6. Either:
   - Paste it into the PWA under **Settings → Cloud Sync Endpoint** for your own device, **OR**
   - Hardcode it into `DEFAULT_SYNC_CONFIG.endpointUrl` at the top of `app/app.js` **before hosting** so every participant automatically syncs to your Sheet without typing anything!
7. Click **Test Cloud Sync** in the app's Settings tab — you'll see a `CheckIns` / `Trials` tab populate automatically on the first entry!

## Option B: Supabase PostgreSQL Collector (Free Tier)
1. Create a free Supabase project and run [`supabase_schema.sql`](./supabase_schema.sql) in the SQL Editor.
2. In the app's **Settings** tab (or `DEFAULT_SYNC_CONFIG` in `app/app.js`), set:
   - **Cloud Backend Mode:** `Supabase REST`
   - **Endpoint URL:** `https://<your-project>.supabase.co/rest/v1/trpm8_events`
   - **Supabase Anon Key:** your project's public `anon` key (safe in client code because RLS restricts `anon` strictly to `INSERT` only!).

---

## Hosting the Mobile PWA Over HTTPS (Required for Mobile GPS & Install-to-Home-Screen)

### Fastest: Netlify Drop (30 seconds)
1. Visit https://app.netlify.com/drop
2. Drag and drop the entire `app/` folder onto the page.
3. Share the generated `https://...netlify.app` link (or plug it straight into a Prolific study URL!).

### Alternative: GitHub Pages (Free Forever)
1. Push `trpm8-osmo-cold/` to a GitHub repository.
2. Go to **Repository Settings → Pages → Deploy from branch → Folder: `/app` (or copy `app/*` to `/docs`)**.
3. Your installable iOS/Android PWA is live at `https://<username>.github.io/<repo>/`.

## Analyzing Collected Data
- **From Google Sheets:** Export `CheckIns` as CSV and `Trials` as CSV (File → Download → CSV), then run:
  ```bash
  python3 analysis/analyze_ema_export.py --input CheckIns.csv Trials.csv
  ```
- **From local JSON exports:**
  ```bash
  python3 analysis/analyze_ema_export.py --input trpm8_study_P-*.json
  ```
