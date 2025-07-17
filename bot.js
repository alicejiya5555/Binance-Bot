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
  const match = cmd.match(/^\/(\w+)(15m|1h|4h|6h|12h)$/);
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

  return vwapArray[vwapArray.length - 1]; // latest VWAP
}

function getKeltnerChannel(candles, emaPeriod = 20, atrPeriod = 14, multiplier = 2) {
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

// 🔧 EMA Helper Function (required by ADOSC)
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

// 📊 KDJ (9,3,3) calculation
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
  const latestD = D[D.length - 1] || 0;
  const J = 3 * latestK - 2 * latestD;

  return {
    k: latestK.toFixed(2),
    d: latestD.toFixed(2),
    j: J.toFixed(2),
  };
}

// 📈 MOMENTUM (MTM) - 7, 14, 20
function getMTM(candles, period) {
  if (candles.length <= period) return 'N/A';

  const currentClose = candles[candles.length - 1].close;
  const pastClose = candles[candles.length - 1 - period].close;
  const mtm = currentClose - pastClose;
  return mtm.toFixed(2);
}

// 📉 ADOSC (Accumulation/Distribution Oscillator)
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

// 🧭 ULTIMATE OSCILLATOR (7,14,28)
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

// 📊 SuperTrend Indicator (ATR Based Trend-Follower)
function getSuperTrend(candles, period = 10, multiplier = 3) {
  const close = candles.map(c => c.close);
  const high = candles.map(c => c.high);
  const low = candles.map(c => c.low);

  const atr = ti.ATR.calculate({ period, high, low, close });
  if (atr.length === 0) return { trend: 'N/A', value: 'N/A' };

  let superTrend = [];
  let upperBand = [];
  let lowerBand = [];
  let trend = 'neutral';

  for (let i = 0; i < close.length; i++) {
    if (i < period) {
      superTrend.push(0);
      upperBand.push(0);
      lowerBand.push(0);
      continue;
    }

    const atrValue = atr[i - 1];
    const hl2 = (high[i] + low[i]) / 2;
    const currentUpper = hl2 + multiplier * atrValue;
    const currentLower = hl2 - multiplier * atrValue;

    upperBand.push(currentUpper);
    lowerBand.push(currentLower);

    if (i === period) {
      superTrend.push(currentUpper);
      trend = 'downtrend';
    } else {
      const prevClose = close[i - 1];
      const prevSuperTrend = superTrend[i - 1];

      if (prevSuperTrend === upperBand[i - 1]) {
        if (close[i] <= currentUpper) {
          superTrend.push(currentUpper);
          trend = 'downtrend';
        } else {
          superTrend.push(currentLower);
          trend = 'uptrend';
        }
      } else {
        if (close[i] >= currentLower) {
          superTrend.push(currentLower);
          trend = 'uptrend';
        } else {
          superTrend.push(currentUpper);
          trend = 'downtrend';
        }
      }
    }
  }

  return {
    trend,
    value: superTrend.length ? superTrend[superTrend.length - 1].toFixed(2) : 'N/A'
  };
}

// 📊 Traders Dynamic Index (TDI)
function getTDI(candles, rsiPeriod = 13, volatilityBand = 34, rsiPriceLine = 2, tradeSignalLine = 7) {
  const close = candles.map(c => c.close);
  
  // Calculate RSI
  const rsi = ti.RSI.calculate({ period: rsiPeriod, values: close });
  if (rsi.length < volatilityBand) return { tdi: 'N/A', trend: 'N/A' };

  // Calculate Bollinger Bands on RSI
  const bb = ti.BollingerBands.calculate({
    period: volatilityBand,
    values: rsi,
    stdDev: rsiPriceLine
  });

  // Calculate Signal Line (MA of RSI)
  const signalLine = ti.SMA.calculate({
    period: tradeSignalLine,
    values: rsi
  });

  if (!bb.length || !signalLine.length) return { tdi: 'N/A', trend: 'N/A' };

  const latestBB = bb[bb.length - 1];
  const latestRSI = rsi[rsi.length - 1];
  const latestSignal = signalLine[signalLine.length - 1];

  let trend = 'neutral';
  if (latestRSI > latestBB.upper) trend = 'overbought';
  else if (latestRSI < latestBB.lower) trend = 'oversold';
  else if (latestRSI > latestSignal) trend = 'bullish';
  else if (latestRSI < latestSignal) trend = 'bearish';

  return {
    tdi: latestRSI.toFixed(2),
    trend,
    upperBand: latestBB.upper.toFixed(2),
    lowerBand: latestBB.lower.toFixed(2),
    signalLine: latestSignal.toFixed(2)
  };
}

// 📊 Heikin Ashi Candles
function getHeikinAshi(candles) {
  if (candles.length < 2) return { trend: 'N/A', candle: 'N/A' };

  const haCandles = [];
  let prevHa = null;

  for (let i = 0; i < candles.length; i++) {
    const current = candles[i];
    
    if (!prevHa) {
      // First candle is same as regular candle
      const haClose = (current.open + current.high + current.low + current.close) / 4;
      prevHa = {
        open: current.open,
        high: current.high,
        low: current.low,
        close: haClose
      };
      haCandles.push(prevHa);
      continue;
    }

    const haClose = (current.open + current.high + current.low + current.close) / 4;
    const haOpen = (prevHa.open + prevHa.close) / 2;
    const haHigh = Math.max(current.high, haOpen, haClose);
    const haLow = Math.min(current.low, haOpen, haClose);

    const haCandle = {
      open: haOpen,
      high: haHigh,
      low: haLow,
      close: haClose
    };

    haCandles.push(haCandle);
    prevHa = haCandle;
  }

  const latestHa = haCandles[haCandles.length - 1];
  let trend = 'neutral';
  
  if (latestHa.close > latestHa.open) {
    trend = latestHa.close > latestHa.high * 0.95 ? 'strong uptrend' : 'uptrend';
  } else if (latestHa.close < latestHa.open) {
    trend = latestHa.close < latestHa.low * 1.05 ? 'strong downtrend' : 'downtrend';
  }

  return {
    trend,
    candle: latestHa,
    close: latestHa.close.toFixed(2)
  };
}

// 📊 Choppiness Index
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

  let marketCondition = 'trending';
  if (ci > 61.8) marketCondition = 'choppy';
  else if (ci < 38.2) marketCondition = 'strong trend';

  return {
    value: Math.min(100, Math.max(0, ci)).toFixed(2),
    condition: marketCondition
  };
}

// --- Indicator Calculations ---
function calculateIndicators(candles) {
  const close = candles.map(c => c.close);
  const high = candles.map(c => c.high);
  const low = candles.map(c => c.low);
  const volume = candles.map(c => c.volume);
  const ichimoku = getIchimoku(candles);
  
  // Helper to safely get last value or NaN if empty
  const lastValue = (arr) => arr.length ? arr.slice(-1)[0] : NaN;

  const macdRaw = ti.MACD.calculate({
    values: close,
    fastPeriod: 3,
    slowPeriod: 10,
    signalPeriod: 16,
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
  const vwap5 = calcVWAP(candles, 5);

  const roc14 = lastValue(ti.ROC.calculate({
    period: 14,
    values: close
  }));

  // 📉 WILLIAMS %R (14)
  function getWilliamsR(candles) {
    const highs = candles.slice(-14).map(c => parseFloat(c[2]));
    const lows = candles.slice(-14).map(c => parseFloat(c[3]));
    const close = parseFloat(candles[candles.length - 1][4]);

    const highestHigh = Math.max(...highs);
    const lowestLow = Math.min(...lows);

    const williamsR = ((highestHigh - close) / (highestHigh - lowestLow)) * -100;
    return williamsR.toFixed(2);
  }

  // 📉 ICHIMOKU (9, 26, 52)
  function getIchimoku(candles) {
    const high = candles.map(c => c.high);
    const low = candles.map(c => c.low);

    const period9 = 9;
    const period26 = 26;
    const period52 = 52;

    if (candles.length < period52) {
      return { conversionLine: 'n/a', baseLine: 'n/a', leadingSpanA: 'n/a', leadingSpanB: 'n/a' };
    }

    const recentHigh9 = Math.max(...high.slice(-period9));
    const recentLow9 = Math.min(...low.slice(-period9));
    const conversionLine = ((recentHigh9 + recentLow9) / 2).toFixed(2);

    const recentHigh26 = Math.max(...high.slice(-period26));
    const recentLow26 = Math.min(...low.slice(-period26));
    const baseLine = ((recentHigh26 + recentLow26) / 2).toFixed(2);

    const leadingSpanA = (((parseFloat(conversionLine) + parseFloat(baseLine)) / 2)).toFixed(2);

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

  // 📊 KDJ indicator calculation
  const kdj = getKDJ(candles);

  const cci7 = lastValue(ti.CCI.calculate({
    period: 7,
    high,
    low,
    close
  }));

  const cci10 = lastValue(ti.CCI.calculate({
    period: 10,
    high,
    low,
    close
  }));

  const cci20 = lastValue(ti.CCI.calculate({
    period: 20,
    high,
    low,
    close
  }));

  const adosc = getADOSC(candles);
  
  // New indicators calculations
  const superTrend = getSuperTrend(candles);
  const tdi = getTDI(candles);
  const heikinAshi = getHeikinAshi(candles);
  const choppinessIndex = getChoppinessIndex(candles);
  
  return {
    sma5: formatNum(lastValue(ti.SMA.calculate({ period: 5, values: close }))),
    sma13: formatNum(lastValue(ti.SMA.calculate({ period: 13, values: close }))),
    sma21: formatNum(lastValue(ti.SMA.calculate({ period: 21, values: close }))),
    sma50: formatNum(lastValue(ti.SMA.calculate({ period: 50, values: close }))),
    sma100: formatNum(lastValue(ti.SMA.calculate({ period: 100, values: close }))),
    sma200: formatNum(lastValue(ti.SMA.calculate({ period: 200, values: close }))),

    ema5: formatNum(lastValue(ti.EMA.calculate({ period: 5, values: close }))),
    ema13: formatNum(lastValue(ti.EMA.calculate({ period: 13, values: close }))),
    ema21: formatNum(lastValue(ti.EMA.calculate({ period: 21, values: close }))),
    ema50: formatNum(lastValue(ti.EMA.calculate({ period: 50, values: close }))),
    ema100: formatNum(lastValue(ti.EMA.calculate({ period: 100, values: close }))),
    ema200: formatNum(lastValue(ti.EMA.calculate({ period: 200, values: close }))),

    wma5: formatNum(lastValue(ti.WMA.calculate({ period: 5, values: close }))),
    wma13: formatNum(lastValue(ti.WMA.calculate({ period: 13, values: close }))),
    wma21: formatNum(lastValue(ti.WMA.calculate({ period: 21, values: close }))),
    wma50: formatNum(lastValue(ti.WMA.calculate({ period: 50, values: close }))),
    wma100: formatNum(lastValue(ti.WMA.calculate({ period: 100, values: close }))),

    macdValue: formatNum(macd.MACD),
    macdSignal: formatNum(macd.signal),
    macdHistogram: formatNum(macd.histogram),

    bbUpper: formatNum(bb.upper),
    bbMiddle: formatNum(bb.middle),
    bbLower: formatNum(bb.lower),

    rsi5: formatNum(lastValue(ti.RSI.calculate({ period: 5, values: close }))),
    rsi14: formatNum(lastValue(ti.RSI.calculate({ period: 14, values: close }))),
   
    atr14: formatNum(atr),

    mfi14: formatNum(lastValue(ti.MFI.calculate({
      high,
      low,
      close,
      volume,
      period: 14
    }))),

    mfi20: formatNum(lastValue(ti.MFI.calculate({
      high,
      low,
      close,
      volume,
      period: 20
    }))),

    williamsR14: formatNum(lastValue(ti.WilliamsR.calculate({
      period: 14,
      high: high,
      low: low,
      close: close
    }))),

    adx14: formatNum(adx),
    pdi14: formatNum(pdi),
    mdi14: formatNum(mdi),

    stochRsiK: formatNum(stochK),
    stochRsiD: formatNum(stochD),

    vwap1: formatNum(vwap1),
    vwap5: formatNum(vwap5),

    kdjK: kdj.k,
    kdjD: kdj.d,
    kdjJ: kdj.j,

    cci7: formatNum(cci7),
    cci10: formatNum(cci10),
    cci20: formatNum(cci20),

    roc14: formatNum(roc14),
    uo: getUltimateOscillator(candles),

    mtm7: getMTM(candles, 7),
    mtm14: getMTM(candles, 14),
    mtm20: getMTM(candles, 20),

    keltner: getKeltnerChannel(candles),

    adosc: isNaN(adosc) ? "N/A" : adosc,

    ichimokuConversion: ichimoku.conversionLine,
    ichimokuBase: ichimoku.baseLine,
    ichimokuSpanA: ichimoku.leadingSpanA,
    ichimokuSpanB: ichimoku.leadingSpanB,
    
    // New indicators
    superTrendValue: superTrend.value,
    superTrendTrend: superTrend.trend,
    
    tdiValue: tdi.tdi,
    tdiTrend: tdi.trend,
    tdiUpperBand: tdi.upperBand,
    tdiLowerBand: tdi.lowerBand,
    tdiSignalLine: tdi.signalLine,
    
    heikinAshiClose: heikinAshi.close,
    heikinAshiTrend: heikinAshi.trend,
    
    choppinessIndex: choppinessIndex.value,
    choppinessCondition: choppinessIndex.condition
  };
}

// --- Output Message Generator ---
function generateOutput(priceData, indicators, name = "Symbol", tfLabel = "Timeframe") {
  const header = 
`📊 ${name} ${tfLabel} Analysis

💰 Price: $${formatNum(priceData.lastPrice)}
📈 24h High: $${formatNum(priceData.highPrice)}
📉 24h Low: $${formatNum(priceData.lowPrice)}
🔁 Change: $${formatNum(priceData.priceChange)} (${priceData.priceChangePercent}%)
🧮 Volume: ${formatNum(priceData.volume)}
💵 Quote Volume: $${formatNum(priceData.quoteVolume)}
🔓 Open Price: $${formatNum(priceData.openPrice)}
⏰ Close Time: ${new Date(priceData.closeTime).toLocaleString('en-UK')}

`;

  const smaSection = 
`📊 Simple Moving Averages (SMA):
 - SMA 5: $${indicators.sma5}
 - SMA 13: $${indicators.sma13}
 - SMA 21: $${indicators.sma21}
 - SMA 50: $${indicators.sma50}
 - SMA 100: $${indicators.sma100}
 - SMA 200: $${indicators.sma200}

`;

  const emaSection =
`📈 Exponential Moving Averages (EMA):
 - EMA 5: $${indicators.ema5}
 - EMA 13: $${indicators.ema13}
 - EMA 21: $${indicators.ema21}
 - EMA 50: $${indicators.ema50}
 - EMA 100: $${indicators.ema100}
 - EMA 200: $${indicators.ema200}

`;

  const wmaSection =
`⚖️ Weighted Moving Averages (WMA):
 - WMA 5: $${indicators.wma5}
 - WMA 13: $${indicators.wma13}
 - WMA 21: $${indicators.wma21}
 - WMA 50: $${indicators.wma50}
 - WMA 100: $${indicators.wma100}

`;

  const macdSection =
`📉 MACD: 3,10,16
 - MACD: ${indicators.macdValue}
 - Signal: ${indicators.macdSignal}
 - Histogram: ${indicators.macdHistogram}

`;

  const bbSection =
`🎯 Bollinger Bands (20, 2 StdDev):
 - Upper Band: $${indicators.bbUpper}
 - Middle Band: $${indicators.bbMiddle}
 - Lower Band: $${indicators.bbLower}

`;

  const rsiSection =
`⚡ Relative Strength Index (RSI):
 - RSI (5): ${indicators.rsi5}
 - RSI (14): ${indicators.rsi14}

`;

  const atrSection = 
`📏 Average True Range (ATR):
 - ATR (14): ${indicators.atr14}

`;

  const adxSection =
`📊 ADX (Trend Strength):
 - ADX (14): ${indicators.adx14}
 - +DI (14): ${indicators.pdi14}
 - -DI (14): ${indicators.mdi14}

`;

  const stochRsiSection =
`📉 Stochastic RSI (14,14,3,3):
 - %K: ${indicators.stochRsiK}
 - %D: ${indicators.stochRsiD}

`;

  const vwapSection =
`🔹 VWAP:
 - VWAP(1): ${indicators.vwap1}
 - VWAP(5): ${indicators.vwap5}

`;

  const mfiSection = 
`💧 Money Flow Index (MFI):
 - MFI (14): ${indicators.mfi14}
 - MFI (20): ${indicators.mfi20}
`;

const williamsSection =
`📉 Williams %R Indicator:
 - Williams %R (14): ${indicators.williamsR14}%
`;

const kdjSection =
`📊 KDJ (9,3,3):
 - K: ${indicators.kdjK}
 - D: ${indicators.kdjD}
 - J: ${indicators.kdjJ}

`;

const cciSection =
`📘 Commodity Channel Index (CCI):
 - CCI (7): ${indicators.cci7}
 - CCI (10): ${indicators.cci10}
 - CCI (20): ${indicators.cci20}

`;

const rocSection =
`📊 Rate of Change (ROC):
 - ROC (14): ${indicators.roc14}%

`;

const uoSection =
`🧭 Ultimate Oscillator:
 - UO (7,14,28): ${indicators.uo}
`;

const mtmSection =
`📈 Momentum (MTM):
 - MTM (7): ${indicators.mtm7}
 - MTM (14): ${indicators.mtm14}
 - MTM (20): ${indicators.mtm20}
`;

const keltnerSection =
`📐 Keltner Channel (20 EMA, 2 ATR):
 - Upper Band: ${indicators.keltner.upper}
 - Middle EMA: ${indicators.keltner.middle}
 - Lower Band: ${indicators.keltner.lower}
`;

const adsocsection = `
📊 ADOSC: ${indicators.adosc}
`;

const ichimokuSection = 
`📊 Ichimoku Cloud:
 - Conversion Line (9): ${indicators.ichimokuConversion}
 - Base Line (26): ${indicators.ichimokuBase}
 - Leading Span A: ${indicators.ichimokuSpanA}
 - Leading Span B: ${indicators.ichimokuSpanB}
`;

// New indicators sections
const superTrendSection =
`📈 SuperTrend (10,3):
 - Value: ${indicators.superTrendValue}
 - Trend: ${indicators.superTrendTrend}

`;

const tdiSection =
`📊 Traders Dynamic Index (TDI):
 - RSI (13): ${indicators.tdiValue}
 - Trend: ${indicators.tdiTrend}
 - Upper Band: ${indicators.tdiUpperBand}
 - Lower Band: ${indicators.tdiLowerBand}
 - Signal Line (7): ${indicators.tdiSignalLine}

`;

const heikinAshiSection =
`🕯 Heikin Ashi:
 - Close: ${indicators.heikinAshiClose}
 - Trend: ${indicators.heikinAshiTrend}

`;

const choppinessSection =
`🌀 Choppiness Index (14):
 - Value: ${indicators.choppinessIndex}
 - Market Condition: ${indicators.choppinessCondition}

`;

  // Your added custom words here:
  const extraNotes =
`
Calculate Values of all thes Indicatotors and Give me Out Put:
📍 Final Signal Summary
📉 What is the overall trend direction? (Bullish, Bearish, or Sideways, positive,Negative, Neutral)
📊 Provide a detailed breakdown of indicator behaviors — RSI, MACD, EMA, Volume, etc.
🔍 Is there any hidden or classic divergence visible (on RSI, MACD, or OBV)?
📈 Are EMAs aligned in a bullish or bearish structure?
🌡 Present a momentum heatmap — Is momentum rising or fading?
📉 Analyze volume and OBV strength — Do they support the price movement?
🧪 Compare current indicators with historically successful setups
⚠️ Scan for breakout or volatility pressure — Are we in a compression or expansion zone?
🔄 After a breakout, is a retest likely? Should we wait for confirmation?
🌬️ Based on news, Twitter, and volume — what’s the real-time sentiment?
🕯 Identify strong candlestick patterns — Engulfing, Doji, Pin Bar, etc.
🔄 Determine whether a reversal or continuation pattern is forming
🌀 Are there any repeating fractal patterns from past cycles?
🐾 Is this setup potentially a bull trap or bear trap?
🧮 Identify key Fibonacci levels — Is the current price near a retracement or extension zone?
🧭 Is the price nearing any known liquidity pool zones?
🛡 Highlight ideal zones for entry, take profit, and stop-loss
🎯 Based on the setup, is TP1, TP2, or TP3 most likely to be hit?
🔁 After taking profit at TP1 or TP2, suggest re-entry levels for the next move
📉 Recommend profitable buy and sell price ranges for this asset
⏳ Compare signals across multiple timeframes (1H, 4H, Daily) — Is there confluence?
🐋 Detect whale movements vs. retail traders — Based on wallet activity or order book flow
🕰 Suggest optimal entry and exit times (based on UTC+07:00 timezone)
📅 Offer a 3-day or weekly forecast — What’s the expected asset behavior?
📰 Is there any upcoming news or event that could impact the market or this asset?
🧠 Suggest the best strategy type for this setup (Scalp, Swing, Position, or News-Driven)
📢 Offer final trading advice — Mindset, Psychology, and Position Sizing
🔁 Is this setup a reversal or continuation opportunity? How clear is the signal?
);

`;

 return header + smaSection + emaSection + wmaSection + macdSection + rsiSection + stochRsiSection + 
        kdjSection + williamsSection + cciSection + rocSection + mtmSection + uoSection + 
        adxSection + bbSection + keltnerSection + atrSection + adsocsection + mfiSection + 
        vwapSection + ichimokuSection + superTrendSection + tdiSection + heikinAshiSection + 
        choppinessSection + extraNotes;
}

// --- Command Handler ---
bot.on("text", async (ctx) => {
  const parsed = parseCommand(ctx.message.text);
  if (!parsed) return ctx.reply("❌ Invalid format. Try `/eth1h`, `/btc15m`, `/link4h`");

  try {
    const { symbol, interval } = parsed;
    const { priceData, candles } = await getBinanceData(symbol, interval);
    const indicators = calculateIndicators(candles);
    
    // Derive friendly names
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
