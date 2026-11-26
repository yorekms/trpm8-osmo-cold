/**
 * Psychrometric & Biophysical Calculations for TRPM8 Osmolarity / Cold Study.
 * Uses August-Roche-Magnus saturation vapor pressure and Stull (2011) wet-bulb
 * approximation to derive true indoor RH and Vapor Pressure Deficit (VPD) from
 * outdoor weather reports — solving the winter indoor-heating RH inversion trap.
 */

/**
 * Saturation vapor pressure e_s(T) in kPa for temperature T in °C.
 * August-Roche-Magnus equation (accurate to <0.4% over -20..50 °C).
 */
export function saturationVaporPressureKPa(tempC) {
  return 0.61094 * Math.exp((17.625 * tempC) / (tempC + 243.04));
}

/**
 * Dew point temperature (°C) given air temperature T (°C) and relative humidity RH (%).
 */
export function dewPointFromTRH(tempC, rhPct) {
  const rh = Math.max(1, Math.min(100, rhPct));
  const gamma = Math.log(rh / 100) + (17.625 * tempC) / (tempC + 243.04);
  return (243.04 * gamma) / (17.625 - gamma);
}

/**
 * Stull (2011) empirical wet-bulb temperature (°C) approximation.
 * Quantifies physical evaporative cooling potential in still air.
 */
export function wetBulbTempC(tempC, rhPct) {
  const rh = Math.max(1, Math.min(100, rhPct));
  return (
    tempC * Math.atan(0.151977 * Math.sqrt(rh + 8.313659)) +
    Math.atan(tempC + rh) -
    Math.atan(rh - 1.676331) +
    0.00391838 * Math.pow(rh, 1.5) * Math.atan(0.023101 * rh) -
    4.686035
  );
}

/**
 * Derive true local microclimate metrics (indoor-corrected when user is indoors).
 *
 * Physics: Unhumidified indoor air retains approximately the outdoor absolute
 * humidity (actual vapor pressure e_a = e_s(T_dew_out)), while warming or cooling
 * to T_in changes e_s(T_in), dramatically altering RH_in and VPD_in.
 *
 * @param {object} params
 * @param {number} params.outdoorTempC
 * @param {number} params.outdoorRhPct
 * @param {number|null} [params.outdoorDewPointC]
 * @param {boolean} params.isIndoors
 * @param {number} [params.indoorTempC=21]
 * @param {boolean} [params.hasHumidifier=false]
 * @param {number} [params.skinTempC=33] - nominal skin/cornea surface temp for surface-to-air VPD
 */
export function computeLocalPsychrometrics({
  outdoorTempC,
  outdoorRhPct,
  outdoorDewPointC = null,
  isIndoors = false,
  indoorTempC = 21.0,
  hasHumidifier = false,
  skinTempC = 33.0,
}) {
  const tdOut =
    outdoorDewPointC !== null && Number.isFinite(outdoorDewPointC)
      ? outdoorDewPointC
      : dewPointFromTRH(outdoorTempC, outdoorRhPct);

  // Actual vapor pressure in outdoor air (kPa)
  const eaOut = saturationVaporPressureKPa(tdOut);
  const esOut = saturationVaporPressureKPa(outdoorTempC);
  const vpdOut = Math.max(0, esOut - eaOut);

  let localTempC = outdoorTempC;
  let localEaKPa = eaOut;

  if (isIndoors) {
    localTempC = Number.isFinite(indoorTempC) ? indoorTempC : 21.0;
    const esIn = saturationVaporPressureKPa(localTempC);
    // Air conditioning condenses moisture if outdoor dew point > indoor temp
    localEaKPa = Math.min(eaOut, esIn * 0.65);
    // Active humidifier raises indoor floor toward ~40% RH if air is very dry
    if (hasHumidifier) {
      localEaKPa = Math.max(localEaKPa, esIn * 0.42);
    }
  }

  const esLocal = saturationVaporPressureKPa(localTempC);
  const localRhPct = Math.max(1, Math.min(100, (localEaKPa / esLocal) * 100));
  const localVpdKPa = Math.max(0, esLocal - localEaKPa);

  // Surface-to-air VPD (driving force for transepidermal water loss & tear evaporation)
  const esSkin = saturationVaporPressureKPa(skinTempC);
  const surfaceVpdKPa = Math.max(0, esSkin - localEaKPa);

  const localWetBulbC = wetBulbTempC(localTempC, localRhPct);
  const evapCoolingPotentialC = Math.max(0, localTempC - localWetBulbC);

  return {
    outdoorTempC: round2(outdoorTempC),
    outdoorRhPct: round1(outdoorRhPct),
    outdoorDewPointC: round2(tdOut),
    outdoorVpdKPa: round3(vpdOut),
    localTempC: round2(localTempC),
    localRhPct: round1(localRhPct),
    localVaporPressureKPa: round3(localEaKPa),
    localVpdKPa: round3(localVpdKPa),
    surfaceVpdKPa: round3(surfaceVpdKPa),
    localWetBulbC: round2(localWetBulbC),
    evapCoolingPotentialC: round2(evapCoolingPotentialC),
  };
}

/**
 * Classify drying stress from local VPD (kPa) for readable user feedback.
 */
export function describeDryingStress(vpdKPa) {
  if (vpdKPa < 0.5) return { label: "Humid / Low Drying", level: "low" };
  if (vpdKPa < 1.2) return { label: "Moderate Drying", level: "moderate" };
  if (vpdKPa < 2.0) return { label: "Dry Air (High TEWL)", level: "high" };
  return { label: "Very Dry / Hyperosmotic Stress", level: "extreme" };
}

function round1(x) {
  return Math.round(x * 10) / 10;
}
function round2(x) {
  return Math.round(x * 100) / 100;
}
function round3(x) {
  return Math.round(x * 1000) / 1000;
}
