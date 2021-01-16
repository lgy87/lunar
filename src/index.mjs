import {
  animals,
  dt_at,
  earthlyBranches,
  heavenlyStems,
  J2000,
  jqmc,
  nutB,
  PI2,
  RAD,
  rmc,
  SSQ_QB,
  SSQ_qiKB,
  SSQ_SB,
  SSQ_suoKB,
  XL0,
  XL0_xzb,
  XL1,
  ymc,
} from "./configs.mjs"

function dt_ext(year, jsd) {
  const dy = (year - 1820) / 100
  return -20 + jsd * dy * dy
}

function dt_calc(y) {
  const y0 = dt_at[dt_at.length - 2]
  const t0 = dt_at[dt_at.length - 1]
  if (y >= y0) {
    const jsd = 31
    if (y > y0 + 100) return dt_ext(y, jsd)
    const v = dt_ext(y, jsd)
    const dv = dt_ext(y0, jsd) - t0
    return v - (dv * (y0 + 100 - y)) / 100
  }

  let i
  for (i = 0; i < dt_at.length; i += 5) {
    if (y < dt_at[i + 5]) break
  }

  const t1 = ((y - dt_at[i]) / (dt_at[i + 5] - dt_at[i])) * 10
  const t2 = t1 * t1
  const t3 = t2 * t1

  return (
    dt_at[i + 1] + dt_at[i + 2] * t1 + dt_at[i + 3] * t2 + dt_at[i + 4] * t3
  )
}

function dt_T(t) {
  return dt_calc(t / 365.2425 + 2000) / 86400.0
}

function gxc_sunLon(t) {
  const squareT = t * t
  const v = -0.043126 + 628.301955 * t - 0.000002732 * squareT
  const e = 0.016708634 - 0.000042037 * t - 0.0000001267 * squareT
  return (-20.49552 * (1 + e * Math.cos(v))) / RAD
}

function XL0_calc(xt, zn, t, n) {
  t /= 10
  let v = 0
  let tn = 1
  const F = XL0[xt]
  const pn = zn * 6 + 1

  const N0 = F[pn + 1] - F[pn]
  for (let i = 0; i < 6; i++, tn *= t) {
    const n1 = F[pn + i]
    const n2 = F[pn + 1 + i]
    const n0 = n2 - n1
    if (!n0) continue
    let N = n2

    if (n >= 0) {
      N = Math.floor((3 * n * n0) / N0 + 0.5) + n1
      if (i) N += 3
      if (N > n2) N = n2
    }
    let j
    let c
    for (j = n1, c = 0; j < N; j += 3) {
      c += F[j] * Math.cos(F[j + 1] + t * F[j + 2])
    }

    v += c * tn
  }

  v /= F[0]

  if (xt == 0) {
    const t2 = t * t
    const t3 = t2 * t
    if (zn == 0) v += (-0.0728 - 2.7702 * t - 1.1019 * t2 - 0.0996 * t3) / RAD
    if (zn == 1) v += (+0.0 + 0.0004 * t + 0.0004 * t2 - 0.0026 * t3) / RAD
    if (zn == 2) v += (-0.002 + 0.0044 * t + 0.0213 * t2 - 0.025 * t3) / 1000000
    return v
  }

  const dv = XL0_xzb[(xt - 1) * 3 + zn]
  if (zn == 0) v += (-3 * t) / RAD
  if (zn == 2) v += dv / 1000000
  else v += dv / RAD
  return v
}

function XL1_calc(zn, t, n) {
  const ob = XL1[zn]
  let t2 = t * t
  let t3 = t2 * t
  let t4 = t3 * t
  const t5 = t4 * t
  const tx = t - 10
  let v = 0

  if (zn == 0) {
    v +=
      (3.81034409 +
        8399.684730072 * t -
        3.319e-5 * t2 +
        3.11e-8 * t3 -
        2.033e-10 * t4) *
      RAD
    v +=
      5028.792262 * t +
      1.1124406 * t2 +
      0.00007699 * t3 -
      0.000023479 * t4 -
      0.0000000178 * t5
    if (tx > 0) v += -0.866 + 1.43 * tx + 0.054 * tx * tx
  }

  t2 /= 1e4
  t3 /= 1e8
  t4 /= 1e8

  n *= 6

  if (n < 0) n = ob[0].length

  for (let i = 0, tn = 1; i < ob.length; i++, tn *= t) {
    const F = ob[i]
    let N = Math.floor((n * F.length) / ob[0].length + 0.5)
    if (i) N += 6
    if (N >= F.length) N = F.length

    let j
    let c
    for (j = 0, c = 0; j < N; j += 6) {
      c +=
        F[j] *
        Math.cos(
          F[j + 1] +
            t * F[j + 2] +
            t2 * F[j + 3] +
            t3 * F[j + 4] +
            t4 * F[j + 5],
        )
    }

    v += c * tn
  }
  if (zn != 2) v /= RAD
  return v
}

function nutationLon2(t) {
  const t2 = t * t
  let dL = 0
  for (let i = 0; i < nutB.length; i += 5) {
    const a = i === 0 ? -1.742 * t : 0
    dL +=
      (nutB[i + 3] + a) * Math.sin(nutB[i] + nutB[i + 1] * t + nutB[i + 2] * t2)
  }
  return dL / 100 / RAD
}

const gxc_moonLon = -3.4e-6

const XL = {
  E_Lon(t, n) {
    return XL0_calc(0, 0, t, n)
  },
  M_Lon(t, n) {
    return XL1_calc(0, t, n)
  },
  E_v(t) {
    const f = 628.307585 * t
    return (
      628.332 +
      21 * Math.sin(1.527 + f) +
      0.44 * Math.sin(1.48 + f * 2) +
      0.129 * Math.sin(5.82 + f) * t +
      0.00055 * Math.sin(4.21 + f) * t * t
    )
  },
  M_v(t) {
    const v =
      8399.71 - 914 * Math.sin(0.7848 + 8328.691425 * t + 0.0001523 * t * t)
    return (
      v -
      179 * Math.sin(2.543 + 15542.7543 * t) +
      160 * Math.sin(0.1874 + 7214.0629 * t) +
      62 * Math.sin(3.14 + 16657.3828 * t) +
      34 * Math.sin(4.827 + 16866.9323 * t) +
      22 * Math.sin(4.9 + 23871.4457 * t) +
      12 * Math.sin(2.59 + 14914.4523 * t) +
      7 * Math.sin(0.23 + 6585.7609 * t) +
      5 * Math.sin(0.9 + 25195.624 * t) +
      5 * Math.sin(2.32 - 7700.3895 * t) +
      5 * Math.sin(3.88 + 8956.9934 * t) +
      5 * Math.sin(0.49 + 7771.3771 * t)
    )
  },
  MS_aLon(t, Mn, Sn) {
    return (
      this.M_Lon(t, Mn) +
      gxc_moonLon -
      this.E_Lon(t, Sn) -
      gxc_sunLon(t) -
      Math.PI
    )
  },
  S_aLon(t, n) {
    return this.E_Lon(t, n) + nutationLon2(t) + gxc_sunLon(t) + Math.PI
  },
  MS_aLon_t(W) {
    const v = 7771.37714500204
    let t = (W + 1.08472) / v
    t += (W - this.MS_aLon(t, 3, 3)) / v

    const v2 = this.M_v(t) - this.E_v(t)
    t += (W - this.MS_aLon(t, 20, 10)) / v2
    t += (W - this.MS_aLon(t, -1, 60)) / v2
    return t
  },
  S_aLon_t(W) {
    let t = (W - 1.75347 - Math.PI) / 628.3319653318
    const ev = this.E_v(t)
    t += (W - this.S_aLon(t, 10)) / ev
    t += (W - this.S_aLon(t, -1)) / ev
    return t
  },
  MS_aLon_t2(W) {
    let v = 7771.37714500204
    let t = (W + 1.08472) / v
    const t2 = t * t
    t -=
      (-0.00003309 * t2 +
        0.10976 * Math.cos(0.784758 + 8328.6914246 * t + 0.000152292 * t2) +
        0.02224 * Math.cos(0.1874 + 7214.0628654 * t - 0.00021848 * t2) -
        0.03342 * Math.cos(4.669257 + 628.307585 * t)) /
      v
    let L =
      this.M_Lon(t, 20) -
      (4.8950632 +
        628.3319653318 * t +
        0.000005297 * t * t +
        0.0334166 * Math.cos(4.669257 + 628.307585 * t) +
        0.0002061 * Math.cos(2.67823 + 628.307585 * t) * t +
        0.000349 * Math.cos(4.6261 + 1256.61517 * t) -
        20.5 / RAD)
    const v2 =
      7771.38 -
      914 * Math.sin(0.7848 + 8328.691425 * t + 0.0001523 * t * t) -
      179 * Math.sin(2.543 + 15542.7543 * t) -
      160 * Math.sin(0.1874 + 7214.0629 * t)
    t += (W - L) / v2
    return t
  },
  S_aLon_t2(W) {
    const v = 628.3319653318
    let t = (W - 1.75347 - Math.PI) / v
    t -=
      (0.000005297 * t * t +
        0.0334166 * Math.cos(4.669257 + 628.307585 * t) +
        0.0002061 * Math.cos(2.67823 + 628.307585 * t) * t) /
      v
    t +=
      (W -
        this.E_Lon(t, 8) -
        Math.PI +
        (20.5 + 17.2 * Math.sin(2.1824 - 33.75705 * t)) / RAD) /
      v
    return t
  },
}

const xAccurate = (W, key) => {
  const t = XL[key](W) * 36525
  return t - dt_T(t) + 8 / 24
}
const qi_accurate = W => xAccurate(W, "S_aLon_t")
const so_accurate = W => xAccurate(W, "MS_aLon_t")

const SSQ = {
  leap: 0,
  ym: [],
  ZQ: [],
  HS: [],
  dx: [],
  Yn: [],
  x_high(W, key, n) {
    const t = XL[`${key}2`](W) * 36525
    const t2 = t - dt_T(t) + 8 / 24
    const v = ((t2 + 0.5) % 1) * 86400

    if (v >= n && v <= 86400 - n) return t2
    return XL[key](W) * 36525 - dt_T(t2) + 8 / 24
  },
  qi_high(W) {
    return this.x_high(W, "S_aLon_t", 1200)
  },
  so_high(W) {
    return this.x_high(W, "MS_aLon_t", 1800)
  },
  calc: function (jd, qs) {
    jd += 2451545
    const B = qs === "气" ? SSQ_qiKB : SSQ_suoKB
    const pc = qs === "气" ? 7 : 14
    const f1 = B[0] - pc
    const f2 = B[B.length - 1] - pc
    const f3 = 2436935
    if (jd < f1 || jd >= f3) {
      if (qs == "气") {
        return Math.floor(
          this.qi_high(
            (Math.floor(((jd + pc - 2451259) / 365.2422) * 24) * Math.PI) / 12,
          ) + 0.5,
        )
      } else {
        return Math.floor(
          this.so_high(
            Math.floor((jd + pc - 2451551) / 29.5306) * Math.PI * 2,
          ) + 0.5,
        )
      }
    }
    let D
    if (jd >= f1 && jd < f2) {
      let i
      for (i = 0; i < B.length; i += 2) {
        if (jd + pc < B[i + 2]) break
      }

      D = B[i] + B[i + 1] * Math.floor((jd + pc - B[i]) / B[i + 1])
      D = Math.floor(D + 0.5)
      if (D === 1683460) D++
      return D - 2451545
    }
    let n
    if (jd >= f2 && jd < f3) {
      if (qs === "气") {
        D = Math.floor(
          this.qi_low(
            (Math.floor(((jd + pc - 2451259) / 365.2422) * 24) * Math.PI) / 12,
          ) + 0.5,
        )
        n = SSQ_QB.substr(Math.floor(((jd - f2) / 365.2422) * 24), 1)
      } else {
        D = Math.floor(
          this.so_low(Math.floor((jd + pc - 2451551) / 29.5306) * Math.PI * 2) +
            0.5,
        )
        n = SSQ_SB.substr(Math.floor((jd - f2) / 29.5306), 1)
      }
      if (n === "1") return D + 1
      if (n === "2") return D - 1
      return D
    }
  },
  calcY: function (jd) {
    const A = this.ZQ
    const B = this.HS

    let W = Math.floor((jd - 355 + 183) / 365.2422) * 365.2422 + 355
    if (this.calc(W, "气") > jd) {
      W -= 365.2422
    }
    let i
    for (i = 0; i < 25; i++) {
      A[i] = this.calc(W + 15.2184 * i, "气")
    }

    A.pe1 = this.calc(W - 15.2, "气")
    A.pe2 = this.calc(W - 30.4, "气")
    let w = this.calc(A[0], "朔")
    if (w > A[0]) w -= 29.53

    for (i = 0; i < 15; i++) {
      B[i] = this.calc(w + 29.5306 * i, "朔")
    }

    this.leap = 0
    for (i = 0; i < 14; i++) {
      this.dx[i] = this.HS[i + 1] - this.HS[i]
      this.ym[i] = i
    }

    const YY = Math.floor((this.ZQ[0] + 10 + 180) / 365.2422) + 2000
    if (YY >= -721 && YY <= -104) {
      const ns = []

      for (i = 0; i < 3; i++) {
        const yy = YY + i - 1

        if (yy >= -721) {
          ns[i] = this.calc(
            1457698 -
              J2000 +
              Math.floor(0.342 + (yy + 721) * 12.368422) * 29.5306,
            "朔",
          )
          ns[i + 3] = "十三"
          ns[i + 6] = 2
        }
        if (yy >= -479) {
          ns[i] = this.calc(
            1546083 -
              J2000 +
              Math.floor(0.5 + (yy + 479) * 12.368422) * 29.5306,
            "朔",
          )
          ns[i + 3] = "十三"
          ns[i + 6] = 2
        }
        if (yy >= -220) {
          ns[i] = this.calc(
            1640641 - J2000 + Math.floor(0.866 + (yy + 220) * 12.369) * 29.5306,
            "朔",
          )
          ns[i + 3] = "后九"
          ns[i + 6] = 11
        }
      }
      for (i = 0; i < 14; i++) {
        let nn
        for (nn = 2; nn >= 0; nn--) {
          if (this.HS[i] >= ns[nn]) break
        }
        const f1 = Math.floor((this.HS[i] - ns[nn] + 15) / 29.5306)

        if (f1 < 12) this.ym[i] = ymc[(f1 + ns[nn + 6]) % 12]
        else this.ym[i] = ns[nn + 3]
      }
      return
    }
    if (B[13] <= A[24]) {
      for (i = 1; B[i + 1] > A[2 * i] && i < 13; i++);
      this.leap = i
      for (; i < 14; i++) this.ym[i]--
    }
    for (i = 0; i < 14; i++) {
      const Dm = this.HS[i] + J2000
      const v2 = this.ym[i]

      let mc = ymc[v2 % 12]
      if (Dm >= 1724360 && Dm <= 1729794) mc = ymc[(v2 + 1) % 12]
      else if (Dm >= 1807724 && Dm <= 1808699) mc = ymc[(v2 + 1) % 12]
      else if (Dm >= 1999349 && Dm <= 1999467) mc = ymc[(v2 + 2) % 12]
      else if (Dm >= 1973067 && Dm <= 1977052) {
        if (v2 % 12 == 0) mc = "正"
        if (v2 == 2) mc = "一"
      }
      if (Dm == 1729794 || Dm == 1808699) mc = "拾贰"
      this.ym[i] = mc
    }
  },
}

const JD = {
  Y: 2000,
  M: 1,
  D: 1,
  h: 12,
  m: 0,
  s: 0,
  JD: function (y, m, d) {
    let n = 0
    let G = 0
    if (y * 372 + m * 31 + Math.floor(d) >= 588829) G = 1
    if (m <= 2) {
      m += 12
      y--
    }
    if (G) {
      n = Math.floor(y / 100)
      n = 2 - n + Math.floor(n / 4)
    }
    return (
      Math.floor(365.25 * (y + 4716)) +
      Math.floor(30.6001 * (m + 1)) +
      d +
      n -
      1524.5
    )
  },
  DD: function (jd) {
    let D = Math.floor(jd + 0.5)
    let F = jd + 0.5 - D
    if (D >= 2299161) {
      const c = Math.floor((D - 1867216.25) / 36524.25)
      D += 1 + c - Math.floor(c / 4)
    }
    D += 1524
    const r = {}
    r.Y = Math.floor((D - 122.1) / 365.25)
    D -= Math.floor(365.25 * r.Y)
    r.M = Math.floor(D / 30.601)
    D -= Math.floor(30.601 * r.M)
    r.D = D
    if (r.M > 13) {
      r.M -= 13
      r.Y -= 4715
    } else {
      r.M -= 1
      r.Y -= 4716
    }
    F *= 24
    r.h = Math.floor(F)
    F -= r.h
    F *= 60
    r.m = Math.floor(F)
    F -= r.m
    F *= 60
    r.s = F
    return r
  },
  toJD: function () {
    return this.JD(
      this.Y,
      this.M,
      this.D + ((this.s / 60 + this.m) / 60 + this.h) / 24,
    )
  },
  setFromJD: function (jd) {
    const r = this.DD(jd)
    this.Y = r.Y
    this.M = r.M
    this.D = r.D
    this.m = r.m
    this.h = r.h
    this.s = r.s
  },
}

function Lunar() {
  this.days = []
  for (let i = 0; i < 31; i++) {
    this.days[i] = {}
  }
}

Lunar.prototype.calc = function (year, month) {
  JD.Y = year
  JD.M = month
  JD.D = 1
  JD.h = 12
  JD.m = 0
  JD.s = 0.1
  const Bd0 = Math.floor(JD.toJD()) - J2000
  JD.M++
  if (JD.M > 12) {
    JD.Y++
    JD.M = 1
  }
  this.daysInMonth = Math.floor(JD.toJD()) - J2000 - Bd0
  this.firstDayInWeek = (Bd0 + J2000 + 1 + 7000000) % 7
  this.year = year
  this.month = month
  const c = year - 1984 + 12000
  this.heavenlyStemsEarthlyBranchesYear =
    heavenlyStems[c % 10] + earthlyBranches[c % 12]
  this.animal = animals[c % 12]
  let D
  for (let i = 0, j = 0; i < this.daysInMonth; i++) {
    const theDay = this.days[i]
    const d0 = Bd0 + i
    theDay.week = (this.firstDayInWeek + i) % 7
    theDay.weekIndex = Math.floor((this.firstDayInWeek + i) / 7)

    JD.setFromJD(d0 + J2000)
    theDay.day = i + 1
    if (!SSQ.ZQ.length || d0 < SSQ.ZQ[0] || d0 >= SSQ.ZQ[24]) SSQ.calcY(d0)
    let mk = Math.floor((d0 - SSQ.HS[0]) / 30)
    if (mk < 13 && SSQ.HS[mk + 1] <= d0) mk++
    theDay.lunarDayIndex = d0 - SSQ.HS[mk]
    theDay.lunarDay = rmc[theDay.lunarDayIndex]
    if (d0 == SSQ.HS[mk] || d0 == Bd0) {
      theDay.lunarMonth = SSQ.ym[mk]
      theDay.lunarDaysInMonth = SSQ.dx[mk]
      theDay.lunarLeap = !!(SSQ.leap && SSQ.leap == mk)
    } else {
      const ob2 = this.days[i - 1]
      theDay.lunarMonth = ob2.lunarMonth
      theDay.lunarDaysInMonth = ob2.lunarDaysInMonth
      theDay.lunarLeap = ob2.lunarLeap
    }
    let qk = Math.floor((d0 - SSQ.ZQ[0] - 7) / 15.2184)
    if (qk < 23 && d0 >= SSQ.ZQ[qk + 1]) qk++
    D = SSQ.ZQ[3] + (d0 < SSQ.ZQ[3] ? -365 : 0) + 365.25 * 16 - 35
    const lyear = Math.floor(D / 365.2422 + 0.5)
    D = SSQ.HS[2]
    for (j = 0; j < 14; j++) {
      if (SSQ.ym[j] != "正" || (SSQ.leap == j && j)) continue
      D = SSQ.HS[j]
      if (d0 < D) {
        D -= 365
        break
      }
    }
    D = D + 5810
    const x = Math.floor(D / 365.2422 + 0.5)
    D = lyear + 12000
    D = x + 12000
    mk = Math.floor((d0 - SSQ.ZQ[0]) / 30.43685)
    if (mk < 12 && d0 >= SSQ.ZQ[2 * mk + 1]) mk++
    D = mk + Math.floor((SSQ.ZQ[12] + 390) / 365.2422) * 12 + 900000
    D = d0 - 6 + 9000000
    mk = Math.floor((d0 - SSQ.ZQ[0] - 15) / 30.43685)
    if (mk < 11 && d0 >= SSQ.ZQ[2 * mk + 2]) mk++
  }

  const jd2 = Bd0 + dt_T(Bd0) - 8 / 24
  let w =
    (Math.floor(((XL.S_aLon(jd2 / 36525, 3) - 0.13) / PI2) * 24) * PI2) / 24
  do {
    const d = qi_accurate(w)
    D = Math.floor(d + 0.5)
    const xn = Math.floor((w / PI2) * 24 + 24000006.01) % 24
    w += PI2 / 24
    if (D >= Bd0 + this.daysInMonth) break
    if (D < Bd0) continue
    this.days[D - Bd0].solarTerm = jqmc[xn]
  } while (D + 12 < Bd0 + this.daysInMonth)
}

export default function lunar(year, month) {
  const lunar = new Lunar()
  lunar.calc(year, month)
  return lunar
}
