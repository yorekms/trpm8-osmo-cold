/**
 * TRPM8 Cold & Hydration Study — Zero-Cost Google Sheets Collector
 *
 * SETUP (3 minutes, $0):
 * 1. Create a new Google Sheet (e.g., "TRPM8 Study Data").
 * 2. Open Extensions -> Apps Script, delete the starter code, and paste this entire file.
 * 3. Click Deploy -> New deployment -> Web app:
 *      - Execute as: Me
 *      - Who has access: Anyone
 * 4. Copy the generated Web App URL (ends in /exec) and paste it into
 *    the TRPM8 PWA -> Settings -> Cloud Sync Webhook URL (or hardcode DEFAULT_WEBHOOK_URL in app.js).
 *
 * PRIVACY: Only anonymized participant IDs, psychrometric values, and sensory scores
 * are accepted. Raw GPS coordinates are never transmitted by the client.
 */

const CHECKIN_HEADERS = [
  "receivedAt", "id", "timestamp", "participantId", "locationMode",
  "indoorTempC", "hasHumidifier", "locationLabel", "outdoorTempC",
  "outdoorRhPct", "outdoorDewPointC", "windKmh", "pressureHpa",
  "localTempC", "localRhPct", "localVpdKPa", "surfaceVpdKPa",
  "evapCoolingPotentialC", "thermalRating", "drySkinRating",
  "dryEyesRating", "wetnessRating", "physicallyWet", "autonomicFlags",
  "waterVolume", "waterRecency", "waterTemp", "waterType",
  "clothing", "activity"
];

const TRIAL_HEADERS = [
  "receivedAt", "id", "timestamp", "participantId", "trialType",
  "localVpdKPa", "localRhPct", "localTempC",
  "mode", "testHand", "bathTempC", "preSkinDeltaC", "forcedChoice",
  "testColdRating", "ctrlColdRating",
  "arm", "timepointMin", "fingerIrTempC", "thermalRating", "shiverUrge", "fingerColdness",
  "occludedSide", "restChoice", "spoonChoice",
  "blinksPer30s", "blinksPerMin", "corneaDryRating", "eyelidCoolRating",
  "notes"
];

function getOrCreateSheet(ss, name, headers) {
  let sh = ss.getSheetByName(name);
  if (!sh) {
    sh = ss.insertSheet(name);
    sh.appendRow(headers);
    sh.setFrozenRows(1);
  }
  return sh;
}

function formatCell(val) {
  if (val === null || val === undefined) return "";
  if (Array.isArray(val)) return val.join(";");
  return val;
}

function doPost(e) {
  try {
    const raw = e && e.postData && e.postData.contents ? e.postData.contents : "{}";
    const body = JSON.parse(raw);
    const recordType = body.recordType || "unknown";
    const rec = body.record || {};

    if (recordType === "ping") {
      return ContentService
        .createTextOutput(JSON.stringify({ ok: true, status: "pong" }))
        .setMimeType(ContentService.MimeType.JSON);
    }

    const ss = SpreadsheetApp.getActiveSpreadsheet();
    const nowIso = new Date().toISOString();

    if (recordType === "checkin") {
      const sh = getOrCreateSheet(ss, "CheckIns", CHECKIN_HEADERS);
      const row = CHECKIN_HEADERS.map((h) => (h === "receivedAt" ? nowIso : formatCell(rec[h])));
      sh.appendRow(row);
    } else if (recordType === "trial") {
      const sh = getOrCreateSheet(ss, "Trials", TRIAL_HEADERS);
      const row = TRIAL_HEADERS.map((h) => (h === "receivedAt" ? nowIso : formatCell(rec[h])));
      sh.appendRow(row);
    } else {
      throw new Error("Unsupported recordType: " + recordType);
    }

    return ContentService
      .createTextOutput(JSON.stringify({ ok: true, recordType, id: rec.id || null }))
      .setMimeType(ContentService.MimeType.JSON);
  } catch (err) {
    return ContentService
      .createTextOutput(JSON.stringify({ ok: false, error: String(err) }))
      .setMimeType(ContentService.MimeType.JSON);
  }
}

function doGet() {
  return ContentService
    .createTextOutput(JSON.stringify({ ok: true, service: "TRPM8 Cold & Hydration Collector" }))
    .setMimeType(ContentService.MimeType.JSON);
}
