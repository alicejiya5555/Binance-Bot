import { Telegraf } from "telegraf";
import axios from "axios";
import ti from "technicalindicators";
import express from "express";

// --- Bot Init ---
const BOT_TOKEN = "7726468556:AAFmVm5S25POmlRXwIRayz1hhbpLP6nDbQ4";
const bot = new Telegraf(BOT_TOKEN);
const PORT = process.env.PORT || 3000;

// --- Utility Functions ---
function parseCommand(command) {
  const cmd = command.toLowerCase();
  const match = cmd.match(/^\/(\w+)(15m|30m|1h|4h|6h|12h)$/);
  if (!match) return null;

  const [, symbolRaw, interval] = match;

  const symbol = symbolRaw === "eth" ? "ETHUSDT"
    : symbolRaw === "btc" ? "BTCUSDT"
    : symbolRaw === "link" ? "LINKUSDT"
    : null;

  if (!symbol) return null;

  return { symbol, interval };
}

function formatNum(num) {
  if (num === undefined || num === null || isNaN(num)) return "N/A";
  return parseFloat(num).toLocaleString("en-US", {
    minimumFractionDigits: 2,
    maximumFractionDigits: 2
  });
}

function calcVWAP(candles, period) {
  let vwapArray = [];
  for (let i = 0; i <= candles.length - period; i++) {
    let slice = candles.slice(i, i + period);
    let cumPV = 0;
    let cumVol = 0;

    for (let bar of slice) {
      const typicalPrice = (parseFloat(bar.high) + parseFloat(bar.low) + parseFloat(bar.close)) / 3;
      const volume = parseFloat(bar.volume);
      cumPV += typicalPrice * volume;
      cumVol += volume;
    }

    vwapArray.push(cumPV / cumVol);
  }

  return vwapArray[vwapArray.length - 1];
}

function getKeltnerChannel(candles, emaPeriod = 20, atrPeriod = 10, multiplier = 2) {
  const close = candles.map(c => c.close);
  const high = candles.map(c => c.high);
  const low = candles.map(c => c.low);

  const emaArray = ti.EMA.calculate({ period: emaPeriod, values: close });
  const atrArray = ti.ATR.calculate({ period: atrPeriod, high, low, close });

  const ema = emaArray.length ? emaArray[emaArray.length - 1] : 0;
  const atr = atrArray.length ? atrArray[atrArray.length - 1] : 0;

  return {
    upper: (ema + multiplier * atr).toFixed(2),
    middle: ema.toFixed(2),
    lower: (ema - multiplier * atr).toFixed(2)
  };
}

// EMA Helper Function
function getEMA(values, period) {
  const k = 2 / (period + 1);
  const emaArray = [];
  let ema = values[0];
  emaArray.push(ema);

  for (let i = 1; i < values.length; i++) {
    ema = values[i] * k + ema * (1 - k);
    emaArray.push(ema);
  }

  return emaArray;
}

// Volume-Weighted MACD (VW-MACD)
function getVWMACD(candles, fastPeriod = 6, slowPeriod = 13, signalPeriod = 5) {
  const close = candles.map(c => c.close);
  const volume = candles.map(c => c.volume);
  
  const vwPrices = close.map((price, i) => price * volume[i]);
  
  const fastEMA = ti.EMA.calculate({ period: fastPeriod, values: vwPrices });
  const slowEMA = ti.EMA.calculate({ period: slowPeriod, values: vwPrices });
  
  const macdLine = [];
  for (let i = 0; i < slowEMA.length; i++) {
    const idx = fastEMA.length - slowEMA.length + i;
    macdLine.push(fastEMA[idx] - slowEMA[i]);
  }
  
  const signalLine = ti.EMA.calculate({ period: signalPeriod, values: macdLine });
  
  const histogram = [];
  for (let i = 0; i < signalLine.length; i++) {
    const idx = macdLine.length - signalLine.length + i;
    histogram.push(macdLine[idx] - signalLine[i]);
  }
  
  return {
    macd: macdLine.length ? macdLine[macdLine.length - 1] : 0,
    signal: signalLine.length ? signalLine[signalLine.length - 1] : 0,
    histogram: histogram.length ? histogram[histogram.length - 1] : 0
  };
}

// Fibonacci Bollinger Bands (FibB)
function getFibonacciBollingerBands(candles, period = 20) {
  const close = candles.map(c => c.close);
  const bb = ti.BollingerBands.calculate({ period, values: close, stdDev: 2 });
  
  if (!bb.length) return { upper: 0, middle: 0, lower: 0 };
  
  const lastBB = bb[bb.length - 1];
  const range = lastBB.upper - lastBB.lower;
  
  return {
    upper: lastBB.upper,
    middle: lastBB.middle,
    lower: lastBB.lower,
    fib0382: (lastBB.middle + range * 0.382).toFixed(2),
    fib0618: (lastBB.middle + range * 0.618).toFixed(2),
    fib1000: lastBB.upper,
    fibNegative0382: (lastBB.middle - range * 0.382).toFixed(2),
    fibNegative0618: (lastBB.middle - range * 0.618).toFixed(2),
    fibNegative1000: lastBB.lower
  };
}

// Relative Volatility Index (RVI)
function getRVI(candles, period = 14) {
  const close = candles.map(c => c.close);
  const high = candles.map(c => c.high);
  const low = candles.map(c => c.low);
  
  const stdevs = [];
  for (let i = period - 1; i < close.length; i++) {
    const slice = close.slice(i - period + 1, i + 1);
    const mean = slice.reduce((a, b) => a + b, 0) / period;
    const variance = slice.reduce((a, b) => a + Math.pow(b - mean, 2), 0) / period;
    stdevs.push(Math.sqrt(variance));
  }
  
  const upChanges = [];
  const downChanges = [];
  
  for (let i = 1; i < stdevs.length; i++) {
    const change = stdevs[i] - stdevs[i - 1];
    if (change > 0) {
      upChanges.push(change);
      downChanges.push(0);
    } else {
      upChanges.push(0);
      downChanges.push(Math.abs(change));
    }
  }
  
  const avgUp = upChanges.reduce((a, b) => a + b, 0) / upChanges.length;
  const avgDown = downChanges.reduce((a, b) => a + b, 0) / downChanges.length;
  
  const rvi = avgDown === 0 ? 100 : 100 - (100 / (1 + (avgUp / avgDown)));
  
  return rvi.toFixed(2);
}

// On-Balance Volume (OBV)
function getOBV(candles) {
  let obv = 0;
  const obvValues = [];
  
  for (let i = 1; i < candles.length; i++) {
    const current = candles[i];
    const previous = candles[i - 1];
    
    if (current.close > previous.close) {
      obv += current.volume;
    } else if (current.close < previous.close) {
      obv -= current.volume;
    }
    
    obvValues.push(obv);
  }
  
  return obvValues.length ? obvValues[obvValues.length - 1] : 0;
}

// Aroon Indicator
function getAroon(candles, period = 14) {
  if (candles.length < period) return { up: 0, down: 0 };
  
  const high = candles.map(c => c.high);
  const low = candles.map(c => c.low);
  
  let aroonUp = 0;
  let aroonDown = 0;
  
  const recentHighs = high.slice(-period);
  const recentLows = low.slice(-period);
  
  const highestHigh = Math.max(...recentHighs);
  const lowestLow = Math.min(...recentLows);
  
  const daysSinceHigh = recentHighs.reverse().findIndex(h => h === highestHigh);
  const daysSinceLow = recentLows.reverse().findIndex(l => l === lowestLow);
  
  aroonUp = ((period - daysSinceHigh) / period) * 100;
  aroonDown = ((period - daysSinceLow) / period) * 100;
  
  return {
    up: aroonUp.toFixed(2),
    down: aroonDown.toFixed(2)
  };
}

// Hull Moving Average (HMA)
function getHMA(candles, period = 9) {
  const close = candles.map(c => c.close);
  
  const halfPeriod = Math.floor(period / 2);
  const wmaHalf = ti.WMA.calculate({ period: halfPeriod, values: close });
  
  const wmaFull = ti.WMA.calculate({ period, values: close });
  
  const rawHMA = [];
  for (let i = 0; i < wmaFull.length; i++) {
    const idx = wmaHalf.length - wmaFull.length + i;
    if (idx >= 0) {
      rawHMA.push(2 * wmaHalf[idx] - wmaFull[i]);
    }
  }
  
  const sqrtPeriod = Math.floor(Math.sqrt(period));
  const hma = ti.WMA.calculate({ period: sqrtPeriod, values: rawHMA });
  
  return hma.length ? hma[hma.length - 1] : 0;
}

// --- Binance Data Fetch ---
async function getBinanceData(symbol, interval) {
  const [priceRes, candlesRes] = await Promise.all([
    axios.get(`https://api.binance.com/api/v3/ticker/24hr?symbol=${symbol}`),
    axios.get(`https://api.binance.com/api/v3/klines?symbol=${symbol}&interval=${interval}&limit=200`)
  ]);

  const priceData = priceRes.data;
  const candles = candlesRes.data.map(c => ({
    time: c[0],
    open: parseFloat(c[1]),
    high: parseFloat(c[2]),
    low: parseFloat(c[3]),
    close: parseFloat(c[4]),
    volume: parseFloat(c[5])
  }));
  
  return { priceData, candles };
}

// KDJ (9,3,3) calculation
function getKDJ(candles) {
  const period = 9;
  const kPeriod = 3;
  const dPeriod = 3;

  const highs = candles.map(c => c.high);
  const lows = candles.map(c => c.low);
  const closes = candles.map(c => c.close);

  const RSV = [];

  for (let i = period - 1; i < closes.length; i++) {
    const highSlice = highs.slice(i - period + 1, i + 1);
    const lowSlice = lows.slice(i - period + 1, i + 1);

    const highestHigh = Math.max(...highSlice);
    const lowestLow = Math.min(...lowSlice);

    const rsv = ((closes[i] - lowestLow) / (highestHigh - lowestLow)) * 100;
    RSV.push(rsv);
  }

  const K = [];
  const D = [];

  K[0] = 50;
  D[0] = 50;

  for (let i = 1; i < RSV.length; i++) {
    K[i] = (2 / 3) * K[i - 1] + (1 / 3) * RSV[i];
    D[i] = (2 / 3) * D[i - 1] + (1 / 3) * K[i];
  }

  const latestK = K[K.length - 1] || 0;
  const latestD = D[K.length - 1] || 0;
  const J = 3 * latestK - 2 * latestD;

  return {
    k: latestK.toFixed(2),
    d: latestD.toFixed(2),
    j: J.toFixed(2),
  };
}

// MOMENTUM (MTM) - 10, 20
function getMTM(candles, period) {
  if (candles.length <= period) return 'N/A';

  const currentClose = candles[candles.length - 1].close;
  const pastClose = candles[candles.length - 1 - period].close;
  const mtm = currentClose - pastClose;
  return mtm.toFixed(2);
}

// ADOSC (Accumulation/Distribution Oscillator)
function getADOSC(candles, fastPeriod = 3, slowPeriod = 10) {
  if (candles.length < slowPeriod) return NaN;

  const adl = [];
  let prevAdl = 0;

  for (let i = 0; i < candles.length; i++) {
    const { high, low, close, volume } = candles[i];
    const hlDiff = high - low;
    const clv = hlDiff === 0 ? 0 : ((close - low) - (high - close)) / hlDiff;
    const moneyFlowVolume = clv * volume;
    const currentAdl = prevAdl + moneyFlowVolume;
    adl.push(currentAdl);
    prevAdl = currentAdl;
  }

  const fastEMA = getEMA(adl, fastPeriod);
  const slowEMA = getEMA(adl, slowPeriod);

  if (!fastEMA.length || !slowEMA.length) return NaN;

  const adosc = fastEMA[fastEMA.length - 1] - slowEMA[slowEMA.length - 1];
  return adosc.toFixed(2);
}

// ULTIMATE OSCILLATOR (7,14,28)
function getUltimateOscillator(candles) {
  if (candles.length < 28) return 'N/A';

  const bp = [];
  const tr = [];

  for (let i = 1; i < candles.length; i++) {
    const curr = candles[i];
    const prev = candles[i - 1];

    const high = curr.high;
    const low = curr.low;
    const closePrev = prev.close;

    const trueLow = Math.min(low, closePrev);
    const trueHigh = Math.max(high, closePrev);

    const buyingPressure = curr.close - trueLow;
    const trueRange = trueHigh - trueLow;

    bp.push(buyingPressure);
    tr.push(trueRange);
  }

  function avg(sumArray, period) {
    const slicedBP = bp.slice(-period);
    const slicedTR = tr.slice(-period);
    const sumBP = slicedBP.reduce((a, b) => a + b, 0);
    const sumTR = slicedTR.reduce((a, b) => a + b, 0);
    return sumTR === 0 ? 0 : sumBP / sumTR;
  }

  const avg7 = avg(bp, 7);
  const avg14 = avg(bp, 14);
  const avg28 = avg(bp, 28);

  const uo = 100 * ((4 * avg7) + (2 * avg14) + avg28) / 7;
  return uo.toFixed(2);
}

// SuperTrend Indicator (ATR Based)
function getSuperTrend(candles, period = 10, multiplier = 3) {
  const close = candles.map(c => parseFloat(c.close));
  const high = candles.map(c => parseFloat(c.high));
  const low = candles.map(c => parseFloat(c.low));

  const atr = ti.ATR.calculate({ period, high, low, close });
  if (atr.length === 0) return { value: 'N/A' };

  let superTrend = [];
  let upperBand = (high[0] + low[0]) / 2;
  let lowerBand = (high[0] + low[0]) / 2;

  for (let i = 0; i < close.length; i++) {
    if (i < period) {
      superTrend.push((high[i] + low[i]) / 2);
      continue;
    }

    const hl2 = (high[i] + low[i]) / 2;
    const currentUpper = hl2 + multiplier * atr[i-1];
    const currentLower = hl2 - multiplier * atr[i-1];

    upperBand = (currentUpper < upperBand || close[i-1] > upperBand) 
      ? currentUpper : upperBand;
    lowerBand = (currentLower > lowerBand || close[i-1] < lowerBand) 
      ? currentLower : lowerBand;

    superTrend.push(
      close[i] > superTrend[i-1] 
        ? Math.max(lowerBand, superTrend[i-1])
        : Math.min(upperBand, superTrend[i-1])
    );
  }

  return {
    value: superTrend.length ? superTrend[superTrend.length - 1].toFixed(2) : 'N/A'
  };
}

// Traders Dynamic Index (TDI)
function getTDI(candles) {
  const close = candles.map(c => c.close);
  const rsi = ti.RSI.calculate({ period: 13, values: close });
  if (rsi.length < 34) return { value: 'N/A', upperBand: 'N/A', lowerBand: 'N/A', signalLine: 'N/A' };

  const bb = ti.BollingerBands.calculate({
    period: 34,
    values: rsi,
    stdDev: 2
  });

  const signalLine = ti.SMA.calculate({
    period: 34,
    values: rsi
  });

  if (!bb.length || !signalLine.length) return { value: 'N/A', upperBand: 'N/A', lowerBand: 'N/A', signalLine: 'N/A' };

  return {
    value: rsi[rsi.length - 1].toFixed(2),
    upperBand: bb[bb.length - 1].upper.toFixed(2),
    lowerBand: bb[bb.length - 1].lower.toFixed(2),
    signalLine: signalLine[signalLine.length - 1].toFixed(2)
  };
}

// Heikin Ashi Candles
function getHeikinAshi(candles) {
  if (candles.length < 2) return { close: 'N/A' };

  const haCandles = [];
  let prevHa = null;

  for (let i = 0; i < candles.length; i++) {
    const current = candles[i];
    
    if (!prevHa) {
      const haClose = (current.open + current.high + current.low + current.close) / 4;
      prevHa = {
        open: current.open,
        close: haClose
      };
      haCandles.push(prevHa);
      continue;
    }

    const haClose = (current.open + current.high + current.low + current.close) / 4;
    const haOpen = (prevHa.open + prevHa.close) / 2;

    const haCandle = {
      open: haOpen,
      close: haClose
    };

    haCandles.push(haCandle);
    prevHa = haCandle;
  }

  return {
    close: haCandles[haCandles.length - 1].close.toFixed(2)
  };
}

// Choppiness Index
function getChoppinessIndex(candles, period = 14) {
  if (candles.length < period + 1) return 'N/A';

  const close = candles.map(c => c.close);
  const high = candles.map(c => c.high);
  const low = candles.map(c => c.low);

  let sumATR = 0;
  let maxHigh = -Infinity;
  let minLow = Infinity;

  for (let i = candles.length - period; i < candles.length; i++) {
    sumATR += Math.max(
      high[i] - low[i],
      Math.abs(high[i] - close[i - 1]),
      Math.abs(low[i] - close[i - 1])
    );

    maxHigh = Math.max(maxHigh, high[i]);
    minLow = Math.min(minLow, low[i]);
  }

  const atrRatio = sumATR / (maxHigh - minLow);
  const ci = 100 * Math.log10(atrRatio) / Math.log10(period);

  return Math.min(100, Math.max(0, ci)).toFixed(2);
}

// Parabolic SAR
function getParabolicSAR(candles, step = 0.02, max = 0.2) {
  const high = candles.map(c => c.high);
  const low = candles.map(c => c.low);
  
  const psar = ti.PSAR.calculate({
    high,
    low,
    step,
    max
  });
  
  return {
    value: psar.length ? psar[psar.length - 1].toFixed(2) : 'N/A'
  };
}

// TRIX Indicator
function getTRIX(candles, period = 10, signalPeriod = 7) {
  const close = candles.map(c => c.close);
  
  const ema1 = ti.EMA.calculate({ period, values: close });
  const ema2 = ti.EMA.calculate({ period, values: ema1 });
  const ema3 = ti.EMA.calculate({ period, values: ema2 });
  
  const trix = [];
  for (let i = 1; i < ema3.length; i++) {
    trix.push((ema3[i] - ema3[i-1]) / ema3[i-1] * 100);
  }
  
  const signal = ti.EMA.calculate({ period: signalPeriod, values: trix });
  
  return {
    value: trix.length ? trix[trix.length - 1].toFixed(4) : 'N/A',
    signal: signal.length ? signal[signal.length - 1].toFixed(4) : 'N/A'
  };
}

// Donchian Channel
function getDonchianChannel(candles, period = 20) {
  if (candles.length < period) return { upper: 'N/A', middle: 'N/A', lower: 'N/A' };
  
  const recentHighs = candles.slice(-period).map(c => c.high);
  const recentLows = candles.slice(-period).map(c => c.low);
  
  const upper = Math.max(...recentHighs);
  const lower = Math.min(...recentLows);
  const middle = (upper + lower) / 2;
  
  return {
    upper: upper.toFixed(2),
    middle: middle.toFixed(2),
    lower: lower.toFixed(2)
  };
}

// Fear & Greed Index (mock implementation)
async function getFearGreedIndex() {
  try {
    return {
      value: Math.floor(Math.random() * 100) + 1,
      classification: ['Extreme Fear', 'Fear', 'Neutral', 'Greed', 'Extreme Greed'][
        Math.floor(Math.random() * 5)
      ]
    };
  } catch (error) {
    return {
      value: 'N/A',
      classification: 'N/A'
    };
  }
}

// ICHIMOKU (9, 26, 52)
function getIchimoku(candles) {
  const high = candles.map(c => c.high);
  const low = candles.map(c => c.low);

  const period9 = 9;
  const period26 = 26;
  const period52 = 52;

  if (candles.length < period52) {
    return { conversionLine: 'N/A', baseLine: 'N/A', leadingSpanA: 'N/A', leadingSpanB: 'N/A' };
  }

  const recentHigh9 = Math.max(...high.slice(-period9));
  const recentLow9 = Math.min(...low.slice(-period9));
  const conversionLine = ((recentHigh9 + recentLow9) / 2).toFixed(2);

  const recentHigh26 = Math.max(...high.slice(-period26));
  const recentLow26 = Math.min(...low.slice(-period26));
  const baseLine = ((recentHigh26 + recentLow26) / 2).toFixed(2);

  const leadingSpanA = ((parseFloat(conversionLine) + parseFloat(baseLine)) / 2).toFixed(2);

  const recentHigh52 = Math.max(...high.slice(-period52));
  const recentLow52 = Math.min(...low.slice(-period52));
  const leadingSpanB = ((recentHigh52 + recentLow52) / 2).toFixed(2);

  return {
    conversionLine,
    baseLine,
    leadingSpanA,
    leadingSpanB
  };
}

// --- Indicator Calculations ---
async function calculateIndicators(candles) {
  const close = candles.map(c => c.close);
  const high = candles.map(c => c.high);
  const low = candles.map(c => c.low);
  const volume = candles.map(c => c.volume);
  const ichimoku = getIchimoku(candles);
  
  const lastValue = (arr) => arr.length ? arr.slice(-1)[0] : NaN;

  const macdRaw = ti.MACD.calculate({
    values: close,
    fastPeriod: 6,
    slowPeriod: 13,
    signalPeriod: 5,
    SimpleMAOscillator: false,
    SimpleMASignal: false
  });
  const macd = lastValue(macdRaw) || { MACD: 0, signal: 0, histogram: 0 };

  const bbRaw = ti.BollingerBands.calculate({
    period: 20,
    values: close,
    stdDev: 2
  });
  const bb = lastValue(bbRaw) || { upper: 0, middle: 0, lower: 0 };

  const atrRaw = ti.ATR.calculate({
    period: 14,
    high,
    low,
    close
  });
  const atr = lastValue(atrRaw);

  const adxData = ti.ADX.calculate({
    period: 14,
    close,
    high,
    low
  });

  const adx = lastValue(adxData)?.adx;
  const pdi = lastValue(adxData)?.pdi;
  const mdi = lastValue(adxData)?.mdi;

  const stochRsiData = ti.StochasticRSI.calculate({
    values: close,
    rsiPeriod: 14,
    stochasticPeriod: 14,
    kPeriod: 3,
    dPeriod: 3
  });

  const stochRsi = lastValue(stochRsiData);
  const stochK = stochRsi?.k;
  const stochD = stochRsi?.d;

  const vwap1 = calcVWAP(candles, 1);
  const vwap3 = calcVWAP(candles, 3);
  const vwap4 = calcVWAP(candles, 4);

  const roc14 = lastValue(ti.ROC.calculate({ period: 14, values: close }));
  const roc25 = lastValue(ti.ROC.calculate({ period: 25, values: close }));

  const vwmacd = getVWMACD(candles, 6, 13, 5);
  const fibBB = getFibonacciBollingerBands(candles);
  const rvi14 = getRVI(candles, 14);
  const rvi10 = getRVI(candles, 10);
  const rviSignal = getRVI(candles, 4);
  const obv = getOBV(candles);
  const aroon = getAroon(candles, 14);
  const hma9 = getHMA(candles, 9);
  const hma14 = getHMA(candles, 14);
  const hma21 = getHMA(candles, 21);

  const kdj = getKDJ(candles);

  const cci14 = lastValue(ti.CCI.calculate({ period: 14, high, low, close }));
  const cci20 = lastValue(ti.CCI.calculate({ period: 20, high, low, close }));

  const adosc = getADOSC(candles);
  
  const superTrend7 = getSuperTrend(candles, 7);
  const superTrend10 = getSuperTrend(candles, 10);
  const superTrend14 = getSuperTrend(candles, 14);
  const tdi = getTDI(candles);
  const heikinAshi = getHeikinAshi(candles);
  const choppinessIndex14 = getChoppinessIndex(candles, 14);
  const choppinessIndex21 = getChoppinessIndex(candles, 21);
  const parabolicSAR = getParabolicSAR(candles, 0.02, 0.2);
  const trix10 = getTRIX(candles, 10, 7);
  const trix14 = getTRIX(candles, 14, 9);
  const donchianChannel20 = getDonchianChannel(candles, 20);
  const donchianChannel14 = getDonchianChannel(candles, 14);
  const fearGreedIndex = await getFearGreedIndex();
  
  return {
    sma10: formatNum(lastValue(ti.SMA.calculate({ period: 10, values: close }))),
    sma20: formatNum(lastValue(ti.SMA.calculate({ period: 20, values: close }))),
    sma50: formatNum(lastValue(ti.SMA.calculate({ period: 50, values: close }))),
    sma200: formatNum(lastValue(ti.SMA.calculate({ period: 200, values: close }))),

    ema9: formatNum(lastValue(ti.EMA.calculate({ period: 9, values: close }))),
    ema21: formatNum(lastValue(ti.EMA.calculate({ period: 21, values: close }))),
    ema50: formatNum(lastValue(ti.EMA.calculate({ period: 50, values: close }))),
    ema200: formatNum(lastValue(ti.EMA.calculate({ period: 200, values: close }))),

    wma8: formatNum(lastValue(ti.WMA.calculate({ period: 8, values: close }))),
    wma20: formatNum(lastValue(ti.WMA.calculate({ period: 20, values: close }))),
    wma50: formatNum(lastValue(ti.WMA.calculate({ period: 50, values: close }))),
    wma100: formatNum(lastValue(ti.WMA.calculate({ period: 100, values: close }))),

    macdValue: formatNum(macd.MACD),
    macdSignal: formatNum(macd.signal),
    macdHistogram: formatNum(macd.histogram),

    bbUpper: formatNum(bb.upper),
    bbMiddle: formatNum(bb.middle),
    bbLower: formatNum(bb.lower),

    rsi3: formatNum(lastValue(ti.RSI.calculate({ period: 3, values: close }))),
    rsi10: formatNum(lastValue(ti.RSI.calculate({ period: 10, values: close }))),
    rsi14: formatNum(lastValue(ti.RSI.calculate({ period: 14, values: close }))),

    atr14: formatNum(atr),

    mfi14: formatNum(lastValue(ti.MFI.calculate({ high, low, close, volume, period: 14 }))),

    williamsR12: formatNum(lastValue(ti.WilliamsR.calculate({ period: 12, high, low, close }))),
    williamsR25: formatNum(lastValue(ti.WilliamsR.calculate({ period: 25, high, low, close }))),

    adx14: formatNum(adx),
    pdi14: formatNum(pdi),
    mdi14: formatNum(mdi),

    stochRsiK: formatNum(stochK),
    stochRsiD: formatNum(stochD),

    vwap1: formatNum(vwap1),
    vwap3: formatNum(vwap3),
    vwap4: formatNum(vwap4),

    kdjK: kdj.k,
    kdjD: kdj.d,
    kdjJ: kdj.j,

    cci14: formatNum(cci14),
    cci20: formatNum(cci20),

    roc14: formatNum(roc14),
    roc25: formatNum(roc25),
    uo: getUltimateOscillator(candles),

    mtm10: getMTM(candles, 10),
    mtm20: getMTM(candles, 20),

    keltner: getKeltnerChannel(candles, 20, 10, 2),

    adosc: isNaN(adosc) ? "N/A" : adosc,

    ichimokuConversion: ichimoku.conversionLine,
    ichimokuBase: ichimoku.baseLine,
    ichimokuSpanA: ichimoku.leadingSpanA,
    ichimokuSpanB: ichimoku.leadingSpanB,
    
    superTrend7: superTrend7.value,
    superTrend10: superTrend10.value,
    superTrend14: superTrend14.value,
    tdi: tdi.value,
    tdiUpperBand: tdi.upperBand,
    tdiLowerBand: tdi.lowerBand,
    tdiSignalLine: tdi.signalLine,
    heikinAshi: heikinAshi.close,
    choppinessIndex14: choppinessIndex14,
    choppinessIndex21: choppinessIndex21,
    
    parabolicSAR: parabolicSAR.value,
    trix10Value: trix10.value,
    trix10Signal: trix10.signal,
    trix14Value: trix14.value,
    trix14Signal: trix14.signal,
    donchianUpper20: donchianChannel20.upper,
    donchianMiddle20: donchianChannel20.middle,
    donchianLower20: donchianChannel20.lower,
    donchianUpper14: donchianChannel14.upper,
    donchianMiddle14: donchianChannel14.middle,
    donchianLower14: donchianChannel14.lower,
    
    vwmacdValue: formatNum(vwmacd.macd),
    vwmacdSignal: formatNum(vwmacd.signal),
    vwmacdHistogram: formatNum(vwmacd.histogram),
    
    fibBBUpper: formatNum(fibBB.upper),
    fibBBMiddle: formatNum(fibBB.middle),
    fibBBLower: formatNum(fibBB.lower),
    fibBB0382: formatNum(fibBB.fib0382),
    fibBB0618: formatNum(fibBB.fib0618),
    fibBB1000: formatNum(fibBB.fib1000),
    fibBBNegative0382: formatNum(fibBB.fibNegative0382),
    fibBBNegative0618: formatNum(fibBB.fibNegative0618),
    fibBBNegative1000: formatNum(fibBB.fibNegative1000),
    
    rvi14: formatNum(rvi14),
    rvi10: formatNum(rvi10),
    rviSignal: formatNum(rviSignal),
    
    aroonUp: aroon.up,
    aroonDown: aroon.down,
    
    hma9: formatNum(hma9),
    hma14: formatNum(hma14),
    hma21: formatNum(hma21)
  };
}

// --- Output Message Generator ---
function generateOutput(priceData, indicators, name = "Symbol", tfLabel = "Timeframe") {
  return `
📊 ${name} ${tfLabel} Analysis

1️⃣ Market Overview
💰 Price: $${formatNum(priceData.lastPrice)}
📈 24h High: $${formatNum(priceData.highPrice)}
📉 24h Low: $${formatNum(priceData.lowPrice)}
🔁 Change: $${formatNum(priceData.priceChange)} (${priceData.priceChangePercent}%)
🧮 Volume: ${formatNum(priceData.volume)}
💵 Quote Volume: $${formatNum(priceData.quoteVolume)}
🔓 Open Price: $${formatNum(priceData.openPrice)}
⏰ Close Time: ${new Date(priceData.closeTime).toLocaleString('en-UK')}

2️⃣ Trend Direction:

📊 Simple Moving Averages (SMA):
 - SMA 10: ${indicators.sma10}
 - SMA 20: ${indicators.sma20}
 - SMA 50: ${indicators.sma50}
 - SMA 200: ${indicators.sma200}

📈 Exponential Moving Averages (EMA):
 - EMA 9: ${indicators.ema9}
 - EMA 21: ${indicators.ema21}
 - EMA 50: ${indicators.ema50}
 - EMA 200: ${indicators.ema200}

⚖️ Weighted Moving Averages (WMA):
 - WMA 8: ${indicators.wma8}
 - WMA 20: ${indicators.wma20}
 - WMA 50: ${indicators.wma50}
 - WMA 100: ${indicators.wma100}

📈 Hull Moving Average:
  (HMA 9): ${indicators.hma9}
  (HMA 14): ${indicators.hma14}
  (HMA 21): ${indicators.hma21}

📊 Ichimoku Cloud:
 - Conversion Line (9): ${indicators.ichimokuConversion}
 - Base Line (26): ${indicators.ichimokuBase}
 - Leading Span A: ${indicators.ichimokuSpanA}
 - Leading Span B: ${indicators.ichimokuSpanB}

📈 SuperTrend:
 - Value(7): ${indicators.superTrend7}
 - Value(10): ${indicators.superTrend10}
 - Value(14): ${indicators.superTrend14}

📈 Parabolic SAR:
 - Step AF Value(0.02): ${indicators.parabolicSAR}
 - Max AF Value(0.20): ${indicators.parabolicSAR}

3️⃣ Momentum Strength

📉 MACD: 6,13,5
 - MACD: ${indicators.macdValue}
 - Signal: ${indicators.macdSignal}
 - Histogram: ${indicators.macdHistogram}

📊 Volume-Weighted MACD (VW-MACD):
 - VW-MACD: ${indicators.vwmacdValue}
 - VW-Signal: ${indicators.vwmacdSignal}
 - VW-Histogram: ${indicators.vwmacdHistogram}

⚡ Relative Strength Index (RSI):
 - RSI (3): ${indicators.rsi3}
 - RSI (10): ${indicators.rsi10}
 - RSI (14): ${indicators.rsi14}

📊 Relative Volatility Index (RVI):
 - RVI (14): ${indicators.rvi14}
 - RVI (10): ${indicators.rvi10}
 - Signal Line(4): ${indicators.rviSignal}

📉 Stochastic RSI (14,3,3)(0.8)level):
 - %K: ${indicators.stochRsiK}
 - %D: ${indicators.stochRsiD}

📊 KDJ (9,3,3):
 - K: ${indicators.kdjK}
 - D: ${indicators.kdjD}
 - J: ${indicators.kdjJ}

📉 Williams %R Indicator:
 - Williams %R (12): ${indicators.williamsR12}
 - Williams %R (25): ${indicators.williamsR25}

📘 Commodity Channel Index (CCI):
 - CCI (14): ${indicators.cci14}
 - CCI (20): ${indicators.cci20}

📊 Rate of Change (ROC):
 - ROC (14): ${indicators.roc14}
 - ROC (25): ${indicators.roc25}

📈 Momentum (MTM):
 - MTM (10): ${indicators.mtm10}
 - MTM (20): ${indicators.mtm20}

🧭 Ultimate Oscillator:
 - UO (7,14,28): ${indicators.uo}

📊 ADX (Trend Strength):
 - ADX (14): ${indicators.adx14}
 - +DI (14): ${indicators.pdi14}
 - -DI (14): ${indicators.mdi14}

📊 Traders Dynamic Index (TDI):
 - RSI (13): ${indicators.tdi}
 - Volatility Bands(34): ${indicators.tdiUpperBand} / ${indicators.tdiLowerBand}
 - Trade Signal Line (34): ${indicators.tdiSignalLine}

4️⃣ Volume & Money Flow

📊 On-Balance Volume (OBV):
 - OBV: ${indicators.obv}

📊 ADOSC: ${indicators.adosc}

💧 Money Flow Index (MFI):
 - MFI (14): ${indicators.mfi14}

📊 Aroon Indicator (14):
 - Aroon Up: ${indicators.aroonUp}
 - Aroon Down: ${indicators.aroonDown}

🔹 VWAP:
 - VWAP(1): ${indicators.vwap1}
 - VWAP(3): ${indicators.vwap3}
 - VWAP(4): ${indicators.vwap4}

5️⃣ Volatility & Range:

🎯 Bollinger Bands (20, 2 StdDev):
 - Upper Band: ${indicators.bbUpper}
 - Middle Band: ${indicators.bbMiddle}
 - Lower Band: ${indicators.bbLower}

📊 Fibonacci Bollinger Bands:
 - Upper (1.0): ${indicators.fibBB1000}
 - Fib 0.618: ${indicators.fibBB0618}
 - Fib 0.382: ${indicators.fibBB0382}
 - Middle: ${indicators.fibBBMiddle}
 - Fib -0.382: ${indicators.fibBBNegative0382}
 - Fib -0.618: ${indicators.fibBBNegative0618}
 - Lower (-1.0): ${indicators.fibBBNegative1000}

📐 Keltner Channel (20 EMA, 2 ATR):
 - Upper Band: ${indicators.keltner.upper}
 - Middle EMA: ${indicators.keltner.middle}
 - Lower Band: ${indicators.keltner.lower}

📏 Average True Range (ATR):
 - ATR (14): ${indicators.atr14}

🕯 Heikin Ashi:
 - Close: ${indicators.heikinAshi}

🌀 Choppiness Index:
 - Value (14): ${indicators.choppinessIndex14}
 - Value (21): ${indicators.choppinessIndex21}
 - Upper Band(61.8): N/A
 - Lower Band(38.2): N/A

📊 TRIX:
 - TRIX(10): ${indicators.trix10Value}
 - TRIX(14): ${indicators.trix14Value}
 - Signal EMA(7): ${indicators.trix10Signal}
 - Signal EMA(9): ${indicators.trix14Signal}

📊 Donchian Channel (20)(14):
 - Upper: ${indicators.donchianUpper20} / ${indicators.donchianUpper14}
 - Middle: ${indicators.donchianMiddle20} / ${indicators.donchianMiddle14}
 - Lower: ${indicators.donchianLower20} / ${indicators.donchianLower14}

Calculate All time frames and Values of all thes Indicatotors and Give me Out Put:
📍 Final Signal Summary
📉 What is the overall trend direction? (Bullish, Bearish, or Sideways, positive,Negative, Neutral)
🌡 Present a momentum heatmap — Is momentum rising or fading?
📉 Analyze volume and OBV strength — Do they support the price movement?
🧪 Compare current indicators with historically successful setups
🔄 Is it any breakout? Or Should we wait for confirmation?
🐾 Over all Is this setup potentially a bull trap or bear trap?
🧮 Identify key Fibonacci levels — Is the current price near a retracement or extension zone?
🧭 Is the price nearing any known liquidity pool zones?
📢 Do you think the price is more likely to decline from here, or is there a greater chance it will rise?
🛡 Highlight ideal zones for entry, take profit, and stop-loss
🎯 Based on the setup, is TP1, TP2, or TP3 most likely to be hit?
🔁 After taking profit at TP1 or TP2, suggest re-entry levels for the next move
⏳ Compare signals across multiple timeframes (1H, 4H, Daily) — Is there confluence?
🐋 Detect whale movements vs. retail traders — Based on wallet activity or order book flow
📅 Offer a 3-day or weekly forecast — What’s the expected asset behavior?
📰 Is there any upcoming news or event that could impact the market or this asset?
📢 Offer final trading advice — Mindset, Psychology, and Position Sizing
`;
}

// --- Command Handler ---
bot.on("text", async (ctx) => {
  const parsed = parseCommand(ctx.message.text);
  if (!parsed) return ctx.reply("❌ Invalid format. Try `/eth1h`, `/btc15m`, `/link4h`");

  try {
    const { symbol, interval } = parsed;
    const { priceData, candles } = await getBinanceData(symbol, interval);
    const indicators = await calculateIndicators(candles);
    
    const name = symbol.replace("USDT", "");
    const tfLabel = interval.toUpperCase();
    
    const message = generateOutput(priceData, indicators, name, tfLabel);
    ctx.reply(message);
  } catch (error) {
    console.error(error);
    ctx.reply("⚠️ Error fetching data. Please try again.");
  }
});

// --- Web Server (keep-alive for Render/Heroku) ---
const app = express();
app.get("/", (req, res) => res.send("Bot is running"));
app.listen(PORT, () => {
  console.log(`Server running on port ${PORT}`);
  bot.launch();
});


