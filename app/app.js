import { computeLocalPsychrometrics, describeDryingStress } from "./psychrometrics.js";

const STORAGE_KEY = "trpm8_cold_study_v1";

// Optional: hardcode your deployed collector URL here before hosting on Netlify/GitHub Pages
// so remote participants sync automatically without touching Settings.
const DEFAULT_SYNC_CONFIG = {
  syncMode: "gas", // "gas" (Google Apps Script) or "supabase"
  endpointUrl: "",
  supabaseAnonKey: "",
};

const THERMAL_LABELS = {
  "-4": "Freezing",
  "-3": "Cold",
  "-2": "Cool",
  "-1": "Slightly Cool",
  "0": "Neutral",
  "1": "Slightly Warm",
  "2": "Warm",
  "3": "Hot",
  "4": "Very Hot",
};

function generateAnonId() {
  const chars = "ABCDEFGHJKLMNPQRSTUVWXYZ23456789";
  let out = "P-";
  for (let i = 0; i < 6; i++) {
    out += chars[Math.floor(Math.random() * chars.length)];
  }
  return out;
}

function loadStore() {
  try {
    const raw = localStorage.getItem(STORAGE_KEY);
    if (raw) {
      const parsed = JSON.parse(raw);
      return {
        participantId: parsed.participantId || generateAnonId(),
        tempUnit: parsed.tempUnit === "C" ? "C" : "F",
        baseDrySkin: parsed.baseDrySkin ?? 3,
        baseDryEye: parsed.baseDryEye ?? 2,
        syncMode: parsed.syncMode || DEFAULT_SYNC_CONFIG.syncMode,
        webhookUrl: parsed.webhookUrl || DEFAULT_SYNC_CONFIG.endpointUrl,
        supabaseAnonKey: parsed.supabaseAnonKey || DEFAULT_SYNC_CONFIG.supabaseAnonKey,
        pendingSync: Array.isArray(parsed.pendingSync) ? parsed.pendingSync : [],
        isDemo: Boolean(parsed.isDemo),
        checkins: Array.isArray(parsed.checkins) ? parsed.checkins : [],
        trials: Array.isArray(parsed.trials) ? parsed.trials : [],
      };
    }
  } catch (err) {
    console.warn("Storage parse failed, resetting:", err);
  }
  return {
    participantId: generateAnonId(),
    tempUnit: "F",
    baseDrySkin: 3,
    baseDryEye: 2,
    syncMode: DEFAULT_SYNC_CONFIG.syncMode,
    webhookUrl: DEFAULT_SYNC_CONFIG.endpointUrl,
    supabaseAnonKey: DEFAULT_SYNC_CONFIG.supabaseAnonKey,
    pendingSync: [],
    isDemo: false,
    checkins: [],
    trials: [],
  };
}

let store = loadStore();
function saveStore() {
  localStorage.setItem(STORAGE_KEY, JSON.stringify(store));
}

// Transient live weather state
let liveWeather = null;
let gelSide = null;
let occSide = null;

// Timer handles
let gelTimerInt = null;
let gelSecLeft = 60;
let bolusTimerInt = null;
let bolusSecElapsed = 0;
let blinkTimerInt = null;
let blinkSecLeft = 30;
let blinkClicks = 0;
let blinkRunning = false;

const $ = (id) => document.getElementById(id);

function cToF(c) {
  return (c * 9) / 5 + 32;
}
function fToC(f) {
  return ((f - 32) * 5) / 9;
}
function fmtTemp(c, decimals = 1) {
  if (!Number.isFinite(c)) return "—";
  if (store.tempUnit === "F") {
    return `${cToF(c).toFixed(decimals)}°F (${c.toFixed(decimals)}°C)`;
  }
  return `${c.toFixed(decimals)}°C (${cToF(c).toFixed(decimals)}°F)`;
}
function fmtTempShort(c, decimals = 1) {
  if (!Number.isFinite(c)) return "—";
  return store.tempUnit === "F"
    ? `${cToF(c).toFixed(decimals)}°F`
    : `${c.toFixed(decimals)}°C`;
}
function fmtDeltaTemp(dc, decimals = 1) {
  if (!Number.isFinite(dc)) return "—";
  return store.tempUnit === "F"
    ? `${(dc * 1.8).toFixed(decimals)}°F`
    : `${dc.toFixed(decimals)}°C`;
}
function updateUnitUI(convertInputs = false) {
  const isF = store.tempUnit === "F";
  const btn = $("btnUnitToggle");
  if (btn) btn.textContent = isF ? "°F (tap→°C)" : "°C (tap→°F)";
  document.querySelectorAll(".unit-lbl").forEach((el) => {
    el.textContent = isF ? "°F" : "°C";
  });
  const ftInp = $("bolusFingerTemp");
  if (ftInp) ftInp.placeholder = isF ? "e.g. 90.3" : "e.g. 32.4";

  if (convertInputs) {
    const bt = parseFloat($("gelBathTemp").value);
    if (Number.isFinite(bt)) {
      $("gelBathTemp").value = (isF ? cToF(bt) : fToC(bt)).toFixed(1);
    }
    const sd = parseFloat($("gelSkinDelta").value);
    if (Number.isFinite(sd)) {
      $("gelSkinDelta").value = (isF ? sd * 1.8 : sd / 1.8).toFixed(1);
    }
    const ft = parseFloat($("bolusFingerTemp").value);
    if (Number.isFinite(ft)) {
      $("bolusFingerTemp").value = (isF ? cToF(ft) : fToC(ft)).toFixed(1);
    }
  } else {
    $("gelBathTemp").value = isF ? "71.6" : "22.0";
  }
  updatePsychrometricDisplay();
  renderDashboard();
}

function showToast(msg, isError = false) {
  const el = $("toast");
  el.textContent = msg;
  el.style.color = isError ? "#f87171" : "#34d399";
  el.style.borderColor = isError ? "#f87171" : "#34d399";
  el.classList.add("show");
  clearTimeout(showToast._t);
  showToast._t = setTimeout(() => el.classList.remove("show"), 2600);
}

// ==================== WEATHER (OPEN-METEO, NO API KEY) ====================

async function fetchOpenMeteo(lat, lon, locationLabel) {
  // Coarse rounding on-device (~11 km privacy grid)
  const rLat = Math.round(lat * 10) / 10;
  const rLon = Math.round(lon * 10) / 10;
  $("weatherStatus").textContent = "Fetching Open-Meteo…";

  const url =
    `https://api.open-meteo.com/v1/forecast?latitude=${rLat}&longitude=${rLon}` +
    `&current=temperature_2m,relative_humidity_2m,dew_point_2m,apparent_temperature,wind_speed_10m,surface_pressure`;

  const resp = await fetch(url);
  if (!resp.ok) throw new Error(`Weather HTTP ${resp.status}`);
  const data = await resp.json();
  const c = data.current;
  if (!c) throw new Error("Missing current weather payload");

  liveWeather = {
    fetchedAt: new Date().toISOString(),
    locationLabel: locationLabel || `${rLat.toFixed(1)}, ${rLon.toFixed(1)}`,
    coarseLat: rLat,
    coarseLon: rLon,
    outdoorTempC: c.temperature_2m,
    outdoorRhPct: c.relative_humidity_2m,
    outdoorDewPointC: c.dew_point_2m,
    apparentTempC: c.apparent_temperature,
    windKmh: c.wind_speed_10m,
    pressureHpa: c.surface_pressure,
  };

  $("weatherStatus").textContent = `Updated · ${liveWeather.locationLabel}`;
  updatePsychrometricDisplay();
  showToast(`Weather loaded (${liveWeather.locationLabel})`);
}

async function lookupCityWeather() {
  const q = $("cityQuery").value.trim();
  if (!q) {
    showToast("Enter a city name first", true);
    return;
  }
  try {
    $("weatherStatus").textContent = "Searching city…";
    const gResp = await fetch(
      `https://geocoding-api.open-meteo.com/v1/search?name=${encodeURIComponent(q)}&count=1&language=en&format=json`
    );
    const gData = await gResp.json();
    if (!gData.results || !gData.results.length) {
      $("weatherStatus").textContent = "City not found";
      showToast("City not found — try a nearby larger city", true);
      return;
    }
    const loc = gData.results[0];
    const label = [loc.name, loc.admin1, loc.country_code].filter(Boolean).join(", ");
    await fetchOpenMeteo(loc.latitude, loc.longitude, label);
  } catch (err) {
    $("weatherStatus").textContent = "Lookup failed";
    showToast("Weather lookup failed (check connection)", true);
  }
}

function fetchGpsWeather() {
  if (!navigator.geolocation) {
    showToast("Geolocation unavailable — use City lookup", true);
    return;
  }
  $("weatherStatus").textContent = "Requesting coarse GPS…";
  navigator.geolocation.getCurrentPosition(
    async (pos) => {
      try {
        await fetchOpenMeteo(pos.coords.latitude, pos.coords.longitude, "Coarse GPS");
      } catch (err) {
        $("weatherStatus").textContent = "Weather fetch error";
        showToast("Could not reach Open-Meteo", true);
      }
    },
    () => {
      $("weatherStatus").textContent = "GPS denied — use City box";
      showToast("GPS denied — type any city instead", true);
    },
    { enableHighAccuracy: false, timeout: 8000, maximumAge: 600000 }
  );
}

function getCurrentPsychrometrics() {
  if (!liveWeather) return null;
  const isIndoors = getSegmentedVal("segLocation") === "indoors";
  const indoorTempC = parseFloat($("rngIndoorTemp").value);
  const hasHumidifier = getSegmentedVal("segHumidifier") === "yes";
  return computeLocalPsychrometrics({
    outdoorTempC: liveWeather.outdoorTempC,
    outdoorRhPct: liveWeather.outdoorRhPct,
    outdoorDewPointC: liveWeather.outdoorDewPointC,
    isIndoors,
    indoorTempC,
    hasHumidifier,
  });
}

function updatePsychrometricDisplay() {
  const psych = getCurrentPsychrometrics();
  if (!psych) return;

  $("mOutTRH").textContent = `${fmtTempShort(psych.outdoorTempC, 1)} / ${psych.outdoorRhPct.toFixed(0)}%`;
  $("mDewPt").textContent = `Dew pt: ${fmtTemp(psych.outdoorDewPointC, 1)}`;

  $("mLocalRh").textContent = `${psych.localRhPct.toFixed(0)}%`;
  $("mLocalTemp").textContent = `Local: ${fmtTemp(psych.localTempC, 1)}`;

  const stress = describeDryingStress(psych.localVpdKPa);
  $("mLocalVpd").textContent = `${psych.localVpdKPa.toFixed(2)} kPa`;
  $("mLocalVpd").className = `val stress-${stress.level}`;
  $("mStressLabel").textContent = stress.label;

  $("mEvapPot").textContent = fmtDeltaTemp(psych.evapCoolingPotentialC, 1);
  const windStr =
    liveWeather.windKmh != null
      ? store.tempUnit === "F"
        ? `${(liveWeather.windKmh * 0.621371).toFixed(1)} mph`
        : `${liveWeather.windKmh} km/h`
      : "—";
  $("mWind").textContent = `Wind: ${windStr}`;
}

// ==================== UI HELPERS ====================

function wireSegmented(containerId, onChange) {
  const box = $(containerId);
  if (!box) return;
  box.addEventListener("click", (e) => {
    const btn = e.target.closest(".chip");
    if (!btn) return;
    box.querySelectorAll(".chip").forEach((c) => c.classList.remove("active"));
    btn.classList.add("active");
    if (onChange) onChange(btn.dataset.val || btn.dataset.trial);
  });
}

function getSegmentedVal(containerId) {
  const active = $(containerId)?.querySelector(".chip.active");
  return active ? active.dataset.val : null;
}

function wireMultiSelect(containerId) {
  const box = $(containerId);
  if (!box) return;
  box.addEventListener("click", (e) => {
    const btn = e.target.closest(".chip");
    if (!btn) return;
    btn.classList.toggle("active");
  });
}

function getMultiSelectFlags(containerId) {
  return Array.from($(containerId)?.querySelectorAll(".chip.active") || []).map(
    (c) => c.dataset.flag
  );
}

function formatThermalLabel(v) {
  const num = parseFloat(v);
  const nearest = String(Math.round(num));
  const txt = THERMAL_LABELS[nearest] || "";
  const sign = num > 0 ? "+" : "";
  return `${sign}${num.toFixed(1)} · ${txt}`;
}

async function postSingleRecord(recordType, record) {
  if (!store.webhookUrl) return false;
  if (store.syncMode === "supabase") {
    const resp = await fetch(store.webhookUrl, {
      method: "POST",
      headers: {
        "Content-Type": "application/json",
        apikey: store.supabaseAnonKey,
        Authorization: `Bearer ${store.supabaseAnonKey}`,
        Prefer: "return=minimal",
      },
      body: JSON.stringify({
        event_id: record.id,
        record_type: recordType,
        participant_id: record.participantId || store.participantId,
        recorded_at: record.timestamp || new Date().toISOString(),
        payload: record,
      }),
    });
    if (!resp.ok) throw new Error(`Supabase HTTP ${resp.status}`);
    return true;
  }

  // Google Apps Script Web App webhook (text/plain avoids CORS preflight)
  await fetch(store.webhookUrl, {
    method: "POST",
    mode: "no-cors",
    headers: { "Content-Type": "text/plain;charset=utf-8" },
    body: JSON.stringify({ recordType, record }),
  });
  return true;
}

function updateSyncBadge() {
  const el = $("syncQueueBadge");
  if (el) {
    const n = store.pendingSync.length;
    el.textContent = n === 0 ? "0 pending" : `${n} queued`;
  }
}

async function syncWebhook(recordType, record) {
  if (!store.webhookUrl) return;
  try {
    if (!navigator.onLine) throw new Error("offline");
    await postSingleRecord(recordType, record);
  } catch (err) {
    store.pendingSync.push({ recordType, record });
    saveStore();
    updateSyncBadge();
    console.warn("Queued record for later cloud sync:", err);
  }
}

async function flushPendingQueue(announce = true) {
  if (!store.webhookUrl) {
    if (announce) showToast("Set Cloud Sync Endpoint URL first", true);
    return;
  }
  if (!store.pendingSync.length) {
    if (announce) showToast("Sync queue is already empty");
    return;
  }
  const remaining = [];
  let sent = 0;
  for (const item of store.pendingSync) {
    try {
      await postSingleRecord(item.recordType, item.record);
      sent++;
    } catch {
      remaining.push(item);
    }
  }
  store.pendingSync = remaining;
  saveStore();
  updateSyncBadge();
  if (announce) {
    showToast(
      remaining.length
        ? `Synced ${sent}; ${remaining.length} still pending`
        : `Synced all ${sent} queued record(s)!`,
      remaining.length > 0
    );
  }
}

async function testCloudSync() {
  if (!store.webhookUrl) {
    showToast("Paste your Google Apps Script or Supabase URL first", true);
    return;
  }
  try {
    if (store.syncMode === "gas") {
      await fetch(store.webhookUrl, {
        method: "POST",
        mode: "no-cors",
        headers: { "Content-Type": "text/plain;charset=utf-8" },
        body: JSON.stringify({ recordType: "ping", record: { participantId: store.participantId } }),
      });
      showToast("Ping sent to Google Apps Script! Log a check-in to see your Sheet populate.");
    } else {
      await postSingleRecord("checkin", {
        id: `ping_${Date.now()}`,
        participantId: store.participantId,
        timestamp: new Date().toISOString(),
        testPing: true,
      });
      showToast("Supabase connection verified (row inserted)!");
    }
  } catch (err) {
    showToast(`Cloud test failed: ${err.message}`, true);
  }
}

// ==================== SAVE CHECK-IN ====================

function logCheckin() {
  const psych = getCurrentPsychrometrics();
  const waterVol = getSegmentedVal("segWaterVol");
  const isIndoors = getSegmentedVal("segLocation") === "indoors";

  const rec = {
    id: `chk_${Date.now()}`,
    timestamp: new Date().toISOString(),
    participantId: store.participantId,
    locationMode: isIndoors ? "indoors" : "outside",
    indoorTempC: isIndoors ? parseFloat($("rngIndoorTemp").value) : null,
    hasHumidifier: isIndoors ? getSegmentedVal("segHumidifier") === "yes" : false,
    locationLabel: liveWeather?.locationLabel || "manual/unfetched",
    outdoorTempC: psych ? psych.outdoorTempC : null,
    outdoorTempF: psych ? Math.round(cToF(psych.outdoorTempC) * 10) / 10 : null,
    outdoorRhPct: psych ? psych.outdoorRhPct : null,
    outdoorDewPointC: psych ? psych.outdoorDewPointC : null,
    windKmh: liveWeather?.windKmh ?? null,
    pressureHpa: liveWeather?.pressureHpa ?? null,
    localTempC: psych ? psych.localTempC : (isIndoors ? parseFloat($("rngIndoorTemp").value) : null),
    localTempF: psych
      ? Math.round(cToF(psych.localTempC) * 10) / 10
      : (isIndoors ? Math.round(cToF(parseFloat($("rngIndoorTemp").value)) * 10) / 10 : null),
    localRhPct: psych ? psych.localRhPct : null,
    localVpdKPa: psych ? psych.localVpdKPa : null,
    surfaceVpdKPa: psych ? psych.surfaceVpdKPa : null,
    evapCoolingPotentialC: psych ? psych.evapCoolingPotentialC : null,
    thermalRating: parseFloat($("rngThermal").value),
    drySkinRating: parseInt($("rngDrySkin").value, 10),
    dryEyesRating: parseInt($("rngDryEyes").value, 10),
    wetnessRating: parseInt($("rngWetness").value, 10),
    physicallyWet: getSegmentedVal("segPhysWet"),
    autonomicFlags: getMultiSelectFlags("segAutonomic"),
    waterVolume: waterVol,
    waterRecency: waterVol === "none" ? "none" : getSegmentedVal("segWaterRecency"),
    waterTemp: waterVol === "none" ? "none" : getSegmentedVal("segWaterTemp"),
    waterType: waterVol === "none" ? "none" : getSegmentedVal("segWaterType"),
    clothing: getSegmentedVal("segClothing"),
    activity: getSegmentedVal("segActivity"),
  };

  store.checkins.push(rec);
  saveStore();
  syncWebhook("checkin", rec);
  renderDashboard();
  showToast(
    psych
      ? `Check-in saved (VPD ${psych.localVpdKPa.toFixed(2)} kPa)`
      : "Check-in saved (tip: fetch weather for VPD!)"
  );
}

// ==================== N-OF-1 TRIALS ====================

function randomizeGel() {
  gelSide = Math.random() < 0.5 ? "LEFT" : "RIGHT";
  const other = gelSide === "LEFT" ? "RIGHT" : "LEFT";
  $("gelAssignment").textContent = `Test = ${gelSide} · Control = ${other}`;
}

function saveGelTrial() {
  const psych = getCurrentPsychrometrics();
  const rawBath = parseFloat($("gelBathTemp").value);
  const rawDelta = parseFloat($("gelSkinDelta").value);
  const bathTempC = store.tempUnit === "F" ? Math.round(fToC(rawBath) * 100) / 100 : rawBath;
  const preSkinDeltaC = store.tempUnit === "F" ? Math.round((rawDelta / 1.8) * 100) / 100 : rawDelta;
  const rec = {
    id: `trl_${Date.now()}`,
    timestamp: new Date().toISOString(),
    participantId: store.participantId,
    trialType: "gel_submersion",
    mode: $("gelMode").value,
    testHand: gelSide || "UNASSIGNED",
    bathTempC,
    bathTempF: Math.round(cToF(bathTempC) * 10) / 10,
    preSkinDeltaC,
    forcedChoice: getSegmentedVal("segGelChoice"),
    testColdRating: parseFloat($("rngGelTestCold").value),
    ctrlColdRating: parseFloat($("rngGelCtrlCold").value),
    notes: $("gelNotes").value.trim(),
    localVpdKPa: psych ? psych.localVpdKPa : null,
  };
  store.trials.push(rec);
  saveStore();
  syncWebhook("trial", rec);
  renderDashboard();
  showToast("Gel submersion trial saved!");
}

function randomizeBolus() {
  const arms = ["warm_water", "warm_saline", "cool_saline", "cool_water"];
  const pick = arms[Math.floor(Math.random() * arms.length)];
  $("bolusArm").value = pick;
  const label = $("bolusArm").selectedOptions[0].textContent.split(":")[0];
  $("bolusArmBadge").textContent = `Assigned: ${label}`;
}

function saveBolusPoint() {
  const psych = getCurrentPsychrometrics();
  const ft = $("bolusFingerTemp").value;
  const rawFt = ft !== "" ? parseFloat(ft) : null;
  const fingerIrTempC =
    rawFt !== null
      ? store.tempUnit === "F"
        ? Math.round(fToC(rawFt) * 100) / 100
        : rawFt
      : null;
  const rec = {
    id: `trl_${Date.now()}`,
    timestamp: new Date().toISOString(),
    participantId: store.participantId,
    trialType: "water_bolus_2x2",
    arm: $("bolusArm").value,
    timepointMin: parseInt($("bolusTimepoint").value, 10),
    fingerIrTempC,
    fingerIrTempF: fingerIrTempC !== null ? Math.round(cToF(fingerIrTempC) * 10) / 10 : null,
    thermalRating: parseFloat($("rngBolusCold").value),
    shiverUrge: parseInt($("rngBolusShiver").value, 10),
    fingerColdness: parseInt($("rngBolusFinger").value, 10),
    localVpdKPa: psych ? psych.localVpdKPa : null,
  };
  store.trials.push(rec);
  saveStore();
  syncWebhook("trial", rec);
  renderDashboard();
  showToast(`Bolus t=${rec.timepointMin}m logged!`);
}

function randomizeOcclusion() {
  occSide = Math.random() < 0.5 ? "LEFT" : "RIGHT";
  const other = occSide === "LEFT" ? "RIGHT" : "LEFT";
  $("occAssignment").textContent = `Occluded = ${occSide} · Bare = ${other}`;
}

function saveOccTrial() {
  const psych = getCurrentPsychrometrics();
  const rec = {
    id: `trl_${Date.now()}`,
    timestamp: new Date().toISOString(),
    participantId: store.participantId,
    trialType: "forearm_occlusion",
    occludedSide: occSide || "UNASSIGNED",
    restChoice: getSegmentedVal("segOccRest"),
    spoonChoice: getSegmentedVal("segOccSpoon"),
    localVpdKPa: psych ? psych.localVpdKPa : null,
    localTempC: psych ? psych.localTempC : null,
  };
  store.trials.push(rec);
  saveStore();
  syncWebhook("trial", rec);
  renderDashboard();
  showToast("Occlusion split trial saved!");
}

function saveBlinkTrial() {
  const psych = getCurrentPsychrometrics();
  const rec = {
    id: `trl_${Date.now()}`,
    timestamp: new Date().toISOString(),
    participantId: store.participantId,
    trialType: "cornea_blink",
    blinksPer30s: blinkClicks,
    blinksPerMin: blinkClicks * 2,
    corneaDryRating: parseInt($("rngCorneaDry").value, 10),
    eyelidCoolRating: parseInt($("rngLidCool").value, 10),
    localVpdKPa: psych ? psych.localVpdKPa : null,
    localRhPct: psych ? psych.localRhPct : null,
    localTempC: psych ? psych.localTempC : null,
  };
  store.trials.push(rec);
  saveStore();
  syncWebhook("trial", rec);
  renderDashboard();
  showToast(`Blink trial saved (${rec.blinksPerMin}/min)!`);
}

// ==================== DASHBOARD & SVG SCATTER CHART ====================

function renderDashboard() {
  $("participantBadge").textContent = `ID: ${store.participantId}`;
  $("demoTag").classList.toggle("hidden", !store.isDemo);

  const checkins = store.checkins;
  const trials = store.trials;
  $("kpiCheckins").textContent = String(checkins.length);
  $("kpiTrials").textContent = String(trials.length);

  const vpds = checkins.map((c) => c.localVpdKPa).filter((v) => Number.isFinite(v));
  $("kpiMeanVpd").textContent = vpds.length
    ? `${(vpds.reduce((a, b) => a + b, 0) / vpds.length).toFixed(2)}`
    : "—";

  const recentWater = checkins.filter(
    (c) =>
      c.waterVolume !== "none" &&
      (c.waterRecency === "lt10" || c.waterRecency === "10_30")
  );
  const noWater = checkins.filter((c) => c.waterVolume === "none");
  if (recentWater.length && noWater.length) {
    const mW = recentWater.reduce((a, b) => a + b.thermalRating, 0) / recentWater.length;
    const mN = noWater.reduce((a, b) => a + b.thermalRating, 0) / noWater.length;
    const diff = mW - mN;
    $("kpiWaterDelta").textContent = `${diff >= 0 ? "+" : ""}${diff.toFixed(2)}`;
  } else {
    $("kpiWaterDelta").textContent = "—";
  }

  renderScatterPlot(checkins);
  renderRecentList(checkins, trials);
}

function renderScatterPlot(checkins) {
  const svg = $("chartScatter");
  const pts = checkins.filter(
    (c) => Number.isFinite(c.localTempC) && Number.isFinite(c.thermalRating)
  );
  if (!pts.length) {
    svg.innerHTML = `<text x="280" y="115" text-anchor="middle" fill="#93a7c6" font-size="13">Log check-ins with weather (or load sample data) to plot Thermal vs. VPD</text>`;
    return;
  }

  const W = 560, H = 220, padL = 44, padR = 20, padT = 18, padB = 34;
  const minT = Math.min(...pts.map((p) => p.localTempC), 5) - 2;
  const maxT = Math.max(...pts.map((p) => p.localTempC), 28) + 2;
  const sx = (t) => padL + ((t - minT) / Math.max(1, maxT - minT)) * (W - padL - padR);
  const sy = (r) => padT + ((4 - r) / 8) * (H - padT - padB);

  let html = "";
  // Grid lines
  for (const r of [-4, -2, 0, 2, 4]) {
    const y = sy(r);
    html += `<line x1="${padL}" y1="${y}" x2="${W - padR}" y2="${y}" stroke="${r === 0 ? "#475569" : "#1e293b"}" stroke-width="1" />`;
    html += `<text x="${padL - 8}" y="${y + 4}" text-anchor="end" fill="#93a7c6" font-size="10">${r > 0 ? "+" + r : r}</text>`;
  }
  const tTicks = 5;
  for (let i = 0; i <= tTicks; i++) {
    const t = minT + ((maxT - minT) * i) / tTicks;
    const x = sx(t);
    const lbl = store.tempUnit === "F" ? `${cToF(t).toFixed(0)}°F` : `${t.toFixed(0)}°C`;
    html += `<line x1="${x}" y1="${padT}" x2="${x}" y2="${H - padB}" stroke="#1e293b" stroke-width="1" />`;
    html += `<text x="${x}" y="${H - 12}" text-anchor="middle" fill="#93a7c6" font-size="10">${lbl}</text>`;
  }
  const axisLbl = store.tempUnit === "F" ? "Local Air Temperature (°F)" : "Local Air Temperature (°C)";
  html += `<text x="${W / 2}" y="${H - 1}" text-anchor="middle" fill="#64748b" font-size="10">${axisLbl}</text>`;

  for (const p of pts) {
    const x = sx(p.localTempC);
    const y = sy(p.thermalRating);
    const vpd = p.localVpdKPa ?? 0.8;
    const fill = vpd >= 1.6 ? "#f87171" : vpd >= 1.0 ? "#fbbf24" : "#38bdf8";
    html += `<circle cx="${x.toFixed(1)}" cy="${y.toFixed(1)}" r="5" fill="${fill}" fill-opacity="0.85" stroke="#0b1220" stroke-width="1"><title>T=${fmtTemp(p.localTempC, 1)}, VPD=${vpd} kPa, Cold=${p.thermalRating}</title></circle>`;
  }
  svg.innerHTML = html;
}

function renderRecentList(checkins, trials) {
  const box = $("recentLogs");
  const combined = [
    ...checkins.map((c) => ({ kind: "checkin", ts: c.timestamp, data: c })),
    ...trials.map((t) => ({ kind: "trial", ts: t.timestamp, data: t })),
  ]
    .sort((a, b) => b.ts.localeCompare(a.ts))
    .slice(0, 12);

  if (!combined.length) {
    box.innerHTML = `<p class="subtext">No logs yet. Record a Check-In or load Sample Pilot Data.</p>`;
    return;
  }

  box.innerHTML = combined
    .map((item) => {
      const dt = new Date(item.ts).toLocaleString([], {
        month: "short",
        day: "numeric",
        hour: "2-digit",
        minute: "2-digit",
      });
      if (item.kind === "checkin") {
        const c = item.data;
        const vpdTxt = Number.isFinite(c.localVpdKPa) ? `VPD ${c.localVpdKPa} kPa` : "No weather";
        const waterTxt =
          c.waterVolume === "none" ? "No recent water" : `Water: ${c.waterVolume} (${c.waterRecency}, ${c.waterType})`;
        return `<div class="log-item"><div class="top"><span>🌡️ Check-In (${c.locationMode})</span><span>${dt}</span></div><div>Thermal <strong>${c.thermalRating}</strong> · DrySkin ${c.drySkinRating}/10 · DryEye ${c.dryEyesRating}/10 · ${vpdTxt} · ${waterTxt}</div></div>`;
      }
      const t = item.data;
      let summary = "";
      if (t.trialType === "gel_submersion") {
        summary = `${t.mode} @ ${fmtTempShort(t.bathTempC, 1)} → ${t.forcedChoice} (Test ${t.testColdRating} vs Ctrl ${t.ctrlColdRating})`;
      } else if (t.trialType === "water_bolus_2x2") {
        summary = `${t.arm} @ t=${t.timepointMin}m · Thermal ${t.thermalRating} · Shiver ${t.shiverUrge}/10 · Finger ${t.fingerColdness}/10`;
      } else if (t.trialType === "forearm_occlusion") {
        summary = `Rest: ${t.restChoice} · Metal spoon: ${t.spoonChoice}`;
      } else {
        summary = `${t.blinksPerMin} blinks/min · CorneaDry ${t.corneaDryRating}/10 · LidCool ${t.eyelidCoolRating}/10`;
      }
      return `<div class="log-item"><div class="top"><span>🧪 ${t.trialType}</span><span>${dt}</span></div><div>${summary}</div></div>`;
    })
    .join("");
}

// ==================== EXPORT & DEMO DATA ====================

function downloadFile(filename, content, mime) {
  const blob = new Blob([content], { type: mime });
  const a = document.createElement("a");
  a.href = URL.createObjectURL(blob);
  a.download = filename;
  a.click();
  setTimeout(() => URL.revokeObjectURL(a.href), 1000);
}

function toCsv(rows) {
  if (!rows.length) return "";
  const keys = Array.from(new Set(rows.flatMap((r) => Object.keys(r))));
  const esc = (v) => {
    const s = Array.isArray(v) ? v.join(";") : v === null || v === undefined ? "" : String(v);
    return /[",\n]/.test(s) ? `"${s.replace(/"/g, '""')}"` : s;
  };
  return [keys.join(","), ...rows.map((r) => keys.map((k) => esc(r[k])).join(","))].join("\n");
}

function loadSamplePilotData() {
  const now = Date.now();
  const demoCheckins = [
    { t: 19, outT: 2, outRh: 75, indoors: true, dry: 7, eye: 6, wet: 1, therm: -2.0, wVol: "none", wRec: "none", wTmp: "none", wTyp: "none" },
    { t: 21, outT: 1, outRh: 80, indoors: true, dry: 8, eye: 7, wet: 0, therm: -1.5, wVol: "large", wRec: "10_30", wTmp: "warm", wTyp: "plain" },
    { t: 21, outT: 18, outRh: 72, indoors: true, dry: 2, eye: 2, wet: 1, therm: 0.0, wVol: "none", wRec: "none", wTmp: "none", wTyp: "none" },
    { t: 20, outT: 16, outRh: 85, indoors: true, dry: 2, eye: 1, wet: 2, therm: 0.0, wVol: "medium", wRec: "10_30", wTmp: "warm", wTyp: "electrolyte" },
    { t: 14, outT: 14, outRh: 30, indoors: false, dry: 7, eye: 6, wet: 2, therm: -2.5, wVol: "none", wRec: "none", wTmp: "none", wTyp: "none" },
    { t: 14, outT: 14, outRh: 82, indoors: false, dry: 2, eye: 1, wet: 3, therm: -1.0, wVol: "none", wRec: "none", wTmp: "none", wTyp: "none" },
    { t: 22, outT: -2, outRh: 70, indoors: true, dry: 8, eye: 7, wet: 0, therm: -1.0, wVol: "large", wRec: "lt10", wTmp: "room", wTyp: "plain" },
    { t: 22, outT: 22, outRh: 65, indoors: false, dry: 2, eye: 2, wet: 1, therm: 0.5, wVol: "none", wRec: "none", wTmp: "none", wTyp: "none" },
  ].map((s, idx) => {
    const p = computeLocalPsychrometrics({
      outdoorTempC: s.outT,
      outdoorRhPct: s.outRh,
      isIndoors: s.indoors,
      indoorTempC: s.t,
    });
    return {
      id: `demo_chk_${idx}`,
      timestamp: new Date(now - (8 - idx) * 3600 * 1000).toISOString(),
      participantId: store.participantId,
      locationMode: s.indoors ? "indoors" : "outside",
      indoorTempC: s.indoors ? s.t : null,
      hasHumidifier: false,
      locationLabel: "Sample Pilot",
      ...p,
      windKmh: 8,
      pressureHpa: 1014,
      thermalRating: s.therm,
      drySkinRating: s.dry,
      dryEyesRating: s.eye,
      wetnessRating: s.wet,
      physicallyWet: "no",
      autonomicFlags: s.therm <= -2 ? ["cold_fingers"] : [],
      waterVolume: s.wVol,
      waterRecency: s.wRec,
      waterTemp: s.wTmp,
      waterType: s.wTyp,
      clothing: "moderate",
      activity: "resting",
    };
  });

  const demoTrials = [
    { id: "demo_t1", timestamp: new Date(now - 7200e3).toISOString(), participantId: store.participantId, trialType: "gel_submersion", mode: "precond_dry_vs_hydrated", testHand: "LEFT", bathTempC: 22.0, preSkinDeltaC: 0.0, forcedChoice: "test_colder", testColdRating: 7.0, ctrlColdRating: 5.5, notes: "Sample dry hand felt sharper cold", localVpdKPa: 1.75 },
    { id: "demo_t2", timestamp: new Date(now - 3600e3).toISOString(), participantId: store.participantId, trialType: "gel_submersion", mode: "hyper_mannitol_vs_iso", testHand: "RIGHT", bathTempC: 22.0, preSkinDeltaC: 0.1, forcedChoice: "test_colder", testColdRating: 6.5, ctrlColdRating: 5.0, notes: "Sample mannitol control", localVpdKPa: 1.60 },
  ];

  store.checkins = demoCheckins;
  store.trials = demoTrials;
  store.isDemo = true;
  saveStore();
  renderDashboard();
  showToast("Loaded sample pilot dataset");
}

// ==================== INIT & EVENT WIRING ====================

function init() {
  // Temperature unit toggle (°F <-> °C)
  $("btnUnitToggle").addEventListener("click", () => {
    store.tempUnit = store.tempUnit === "F" ? "C" : "F";
    saveStore();
    updateUnitUI(true);
    showToast(`Display switched to °${store.tempUnit} (exports preserve SI °C)`);
  });

  // Bottom tabs
  document.querySelectorAll(".nav-btn").forEach((btn) => {
    btn.addEventListener("click", () => {
      document.querySelectorAll(".nav-btn").forEach((b) => b.classList.remove("active"));
      document.querySelectorAll(".tab-panel").forEach((p) => p.classList.remove("active"));
      btn.classList.add("active");
      $(btn.dataset.tab).classList.add("active");
    });
  });

  // Weather
  $("btnGpsWeather").addEventListener("click", fetchGpsWeather);
  $("btnCityWeather").addEventListener("click", lookupCityWeather);

  wireSegmented("segLocation", (val) => {
    $("indoorControls").classList.toggle("hidden", val !== "indoors");
    updatePsychrometricDisplay();
  });
  wireSegmented("segHumidifier", () => updatePsychrometricDisplay());
  $("rngIndoorTemp").addEventListener("input", (e) => {
    const c = parseFloat(e.target.value);
    $("valIndoorTemp").textContent = c.toFixed(1);
    $("valIndoorTempF").textContent = Math.round((c * 9) / 5 + 32);
    updatePsychrometricDisplay();
  });

  // Qualia sliders
  $("rngThermal").addEventListener("input", (e) => {
    $("valThermal").textContent = formatThermalLabel(e.target.value);
  });
  $("rngDrySkin").addEventListener("input", (e) => ($("valDrySkin").textContent = `${e.target.value} / 10`));
  $("rngDryEyes").addEventListener("input", (e) => ($("valDryEyes").textContent = `${e.target.value} / 10`));
  $("rngWetness").addEventListener("input", (e) => ($("valWetness").textContent = `${e.target.value} / 10`));
  wireSegmented("segPhysWet");
  wireMultiSelect("segAutonomic");

  // Water controls
  wireSegmented("segWaterVol", (val) => {
    $("waterDetailBox").classList.toggle("hidden", val === "none");
  });
  wireSegmented("segWaterRecency");
  wireSegmented("segWaterTemp");
  wireSegmented("segWaterType");
  wireSegmented("segClothing");
  wireSegmented("segActivity");
  $("btnLogCheckin").addEventListener("click", logCheckin);

  // Trials navigation
  wireSegmented("segTrialSelect", (trial) => {
    document.querySelectorAll(".trial-pane").forEach((p) => p.classList.add("hidden"));
    $(`trial-${trial}`).classList.remove("hidden");
  });

  // Trial A: Gel
  $("btnRandGel").addEventListener("click", randomizeGel);
  wireSegmented("segGelChoice");
  $("rngGelTestCold").addEventListener("input", (e) => ($("valGelTestCold").textContent = `${e.target.value} / 10`));
  $("rngGelCtrlCold").addEventListener("input", (e) => ($("valGelCtrlCold").textContent = `${e.target.value} / 10`));
  $("btnGelTimerStart").addEventListener("click", () => {
    clearInterval(gelTimerInt);
    gelSecLeft = 60;
    $("gelTimerDisplay").textContent = "01:00";
    gelTimerInt = setInterval(() => {
      gelSecLeft--;
      const m = String(Math.floor(gelSecLeft / 60)).padStart(2, "0");
      const s = String(gelSecLeft % 60).padStart(2, "0");
      $("gelTimerDisplay").textContent = `${m}:${s}`;
      if (gelSecLeft <= 0) {
        clearInterval(gelTimerInt);
        showToast("60s gel immersion complete — rate both hands now!");
      }
    }, 1000);
  });
  $("btnGelTimerReset").addEventListener("click", () => {
    clearInterval(gelTimerInt);
    gelSecLeft = 60;
    $("gelTimerDisplay").textContent = "01:00";
  });
  $("btnSaveGel").addEventListener("click", saveGelTrial);

  // Trial B: Bolus
  $("btnRandBolus").addEventListener("click", randomizeBolus);
  $("rngBolusCold").addEventListener("input", (e) => ($("valBolusCold").textContent = e.target.value));
  $("rngBolusShiver").addEventListener("input", (e) => ($("valBolusShiver").textContent = `${e.target.value} / 10`));
  $("rngBolusFinger").addEventListener("input", (e) => ($("valBolusFinger").textContent = `${e.target.value} / 10`));
  $("btnBolusTimerStart").addEventListener("click", () => {
    clearInterval(bolusTimerInt);
    bolusTimerInt = setInterval(() => {
      bolusSecElapsed++;
      const m = String(Math.floor(bolusSecElapsed / 60)).padStart(2, "0");
      const s = String(bolusSecElapsed % 60).padStart(2, "0");
      $("bolusTimerDisplay").textContent = `${m}:${s}`;
      if ([300, 900, 1800].includes(bolusSecElapsed)) {
        showToast(`Bolus checkpoint reached (${bolusSecElapsed / 60} min)!`);
      }
    }, 1000);
  });
  $("btnBolusTimerReset").addEventListener("click", () => {
    clearInterval(bolusTimerInt);
    bolusSecElapsed = 0;
    $("bolusTimerDisplay").textContent = "00:00";
  });
  $("btnSaveBolusPoint").addEventListener("click", saveBolusPoint);

  // Trial C: Occlusion
  $("btnRandOcc").addEventListener("click", randomizeOcclusion);
  wireSegmented("segOccRest");
  wireSegmented("segOccSpoon");
  $("btnSaveOcc").addEventListener("click", saveOccTrial);

  // Trial D: Blink
  $("rngCorneaDry").addEventListener("input", (e) => ($("valCorneaDry").textContent = `${e.target.value} / 10`));
  $("rngLidCool").addEventListener("input", (e) => ($("valLidCool").textContent = `${e.target.value} / 10`));
  $("btnBlinkStart").addEventListener("click", () => {
    clearInterval(blinkTimerInt);
    blinkSecLeft = 30;
    blinkClicks = 0;
    blinkRunning = true;
    $("blinkTime").textContent = "30";
    $("blinkCount").textContent = "0";
    blinkTimerInt = setInterval(() => {
      blinkSecLeft--;
      $("blinkTime").textContent = String(blinkSecLeft);
      if (blinkSecLeft <= 0) {
        clearInterval(blinkTimerInt);
        blinkRunning = false;
        showToast(`30s done: ${blinkClicks} blinks (${blinkClicks * 2}/min)`);
      }
    }, 1000);
  });
  $("btnBlinkReset").addEventListener("click", () => {
    clearInterval(blinkTimerInt);
    blinkSecLeft = 30;
    blinkClicks = 0;
    blinkRunning = false;
    $("blinkTime").textContent = "30";
    $("blinkCount").textContent = "0";
  });
  $("btnBlinkTap").addEventListener("click", () => {
    if (!blinkRunning && blinkSecLeft === 30) {
      $("btnBlinkStart").click();
    }
    if (blinkSecLeft > 0) {
      blinkClicks++;
      $("blinkCount").textContent = String(blinkClicks);
    }
  });
  $("btnSaveBlink").addEventListener("click", saveBlinkTrial);

  // Export & Demo
  $("btnExportJson").addEventListener("click", () => {
    const payload = {
      schemaVersion: "1.0",
      exportedAt: new Date().toISOString(),
      participant: {
        participantId: store.participantId,
        baseDrySkin: store.baseDrySkin,
        baseDryEye: store.baseDryEye,
      },
      checkins: store.checkins,
      trials: store.trials,
    };
    downloadFile(
      `trpm8_study_${store.participantId}_${Date.now()}.json`,
      JSON.stringify(payload, null, 2),
      "application/json"
    );
  });
  $("btnExportCheckinsCsv").addEventListener("click", () => {
    if (!store.checkins.length) return showToast("No check-ins to export", true);
    downloadFile(`trpm8_checkins_${store.participantId}.csv`, toCsv(store.checkins), "text/csv");
  });
  $("btnExportTrialsCsv").addEventListener("click", () => {
    if (!store.trials.length) return showToast("No trials to export", true);
    downloadFile(`trpm8_trials_${store.participantId}.csv`, toCsv(store.trials), "text/csv");
  });
  $("btnLoadDemo").addEventListener("click", loadSamplePilotData);
  $("btnClearData").addEventListener("click", () => {
    if (!confirm("Clear all locally saved check-ins and trials?")) return;
    store.checkins = [];
    store.trials = [];
    store.isDemo = false;
    saveStore();
    renderDashboard();
    showToast("Local data cleared");
  });

  // Settings
  $("inpParticipantId").value = store.participantId;
  $("inpBaseDrySkin").value = store.baseDrySkin;
  $("inpBaseDryEye").value = store.baseDryEye;
  $("selSyncMode").value = store.syncMode;
  $("inpWebhookUrl").value = store.webhookUrl;
  $("inpSupabaseKey").value = store.supabaseAnonKey;
  $("supabaseKeyBox").classList.toggle("hidden", store.syncMode !== "supabase");
  updateSyncBadge();

  $("selSyncMode").addEventListener("change", (e) => {
    $("supabaseKeyBox").classList.toggle("hidden", e.target.value !== "supabase");
  });
  $("btnSaveSettings").addEventListener("click", () => {
    store.participantId = $("inpParticipantId").value.trim() || generateAnonId();
    store.baseDrySkin = parseInt($("inpBaseDrySkin").value, 10) || 0;
    store.baseDryEye = parseInt($("inpBaseDryEye").value, 10) || 0;
    store.syncMode = $("selSyncMode").value;
    store.webhookUrl = $("inpWebhookUrl").value.trim();
    store.supabaseAnonKey = $("inpSupabaseKey").value.trim();
    saveStore();
    renderDashboard();
    updateSyncBadge();
    showToast("Settings saved");
  });
  $("btnTestSync").addEventListener("click", testCloudSync);
  $("btnFlushSync").addEventListener("click", () => flushPendingQueue(true));
  window.addEventListener("online", () => flushPendingQueue(false));

  // Service worker registration for offline PWA install
  if ("serviceWorker" in navigator) {
    navigator.serviceWorker.register("./sw.js").catch(() => {});
  }

  updateUnitUI(false);
}

init();
