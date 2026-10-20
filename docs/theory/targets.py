"""Independent computation of numeric test targets for docs/theory/*.md.
Uses only closed forms (no cvtsim code). Run: python docs/theory/targets.py"""
import numpy as np
from scipy.special import ellipk
from scipy.optimize import brentq

out = {}
# 1. pendulum (L=1 m, g=9.81)
g, L = 9.81, 1.0
for th0 in (10, 60, 120, 170):
    k = np.sin(np.radians(th0) / 2)
    out[f"pend_T_{th0}deg"] = 4*np.sqrt(L/g)*ellipk(k**2)
# torque-free symmetric top I1=I2=1, I3=2, omega_body=(1,0,3)
I1, I3 = 1.0, 2.0
w = np.array([1.0, 0.0, 3.0])
out["top_body_prec_rate"] = (I3-I1)/I1*w[2]
Lmag = np.linalg.norm([I1*w[0], I1*w[1], I3*w[2]])
out["top_space_prec_rate"] = Lmag/I1
out["top_nutation_angle_deg"] = np.degrees(np.arccos(I3*w[2]/Lmag))
# 2. gen-alpha params rho_inf=0.8
r = 0.8
am = (2*r-1)/(r+1); af = r/(r+1); gam = 0.5-am+af; bet = 0.25*(1-am+af)**2
out.update(ga_am=am, ga_af=af, ga_gamma=gam, ga_beta=bet)
# 3/4. belt / contact
beta = np.radians(45.0); mu = 0.3
out["capstan_ratio_pi"] = np.exp(mu/np.sin(beta)*np.pi)
out["Fa_T10_R0p05"] = 10*np.cos(beta)/(2*mu*0.05)
EA = 243.5e3
out["stretch_T500_L1"] = 500*1.0/EA
# 5. CVT
n, mfw = 3, 0.275+0.0387; rw0 = 0.04445; k_p, F0p = 13208.0, 683.04
for kap in (1.0, 0.5):
    wE = np.sqrt(F0p/(n*mfw*rw0*kap)); out[f"engage_rpm_kappa{kap}"] = wE*60/(2*np.pi)
    wN = np.sqrt(k_p/(n*mfw*kap**2)); out[f"negstiff_rpm_kappa{kap}"] = wN*60/(2*np.pi)
kt, th0, rc, phi = 3.930, np.radians(30), 0.0381, np.radians(45)
out["sec_torsion_axial_preload"] = kt*th0/(rc*np.tan(phi))
out["sec_linear_preload_derived"] = 435.63 - out["sec_torsion_axial_preload"]
out["sec_thrust_T20_x0"] = 435.63 + (20/2)/(rc*np.tan(phi))
# geometry: measured hardware (05 sec 2), belt-length closure
def Lbelt(r1, r2, C):
    ph = np.arcsin((r2-r1)/C)
    return 2*C*np.cos(ph) + np.pi*(r1+r2) + 2*ph*(r2-r1)
Lb, Cc, b12 = 0.88009, 0.232, np.radians(12.0)
rs_of = lambda r: brentq(lambda q: Lbelt(r, q, Cc) - Lb, 1e-6, r + 0.999*Cc)
rp_lo = 0.0315; rs_lo = rs_of(rp_lo)
rp_hi = brentq(lambda r: rs_of(r)/r - 0.8, rp_lo, 0.0914)
# secondary-limited case: real belt BD52-2179-S (895 mm), hub 26.416 mm, secondary max 97.03 mm
Lb2 = 0.895
rs_of2 = lambda r: brentq(lambda q: Lbelt(r, q, Cc) - Lb2, 1e-6, r + 0.999*Cc)
rp_lo2 = brentq(lambda r: rs_of2(r) - 0.09703, 0.026416, 0.0914)
rp_hi2 = brentq(lambda r: rs_of2(r)/r - 0.8, rp_lo2, 0.0914)
out.update(geo2_rp_low=rp_lo2, geo2_ratio_low=0.09703/rp_lo2, geo2_rp_max=rp_hi2,
           geo2_rs_min=rs_of2(rp_hi2), geo2_x_p_stop=2*np.tan(b12)*(rp_hi2 - rp_lo2))
out.update(geo_rs_low=rs_lo, geo_ratio_low=rs_lo/rp_lo, geo_rp_max=rp_hi, geo_rs_min=rs_of(rp_hi),
           geo_x_p_stop=2*np.tan(b12)*(rp_hi - rp_lo), geo_x_s0=2*np.tan(b12)*(0.09703 - rs_lo),
           geo_rp_at_ratio3=brentq(lambda r: rs_of(r)/r - 3.0, 0.02, 0.05))
for k_, v in out.items():
    print(f"{k_:32s} {v:.6g}")

# --- drivetrain steady state: T(w)=a-b*w, ratio R, G, rt, Crr, m, c_a
a_, b_, R_, G_, rt_, Crr_, m_, ca_ = 40.0, 0.05, 1.0, 8.1, 0.2921, 0.015, 250.0, 0.5
k_ = R_*G_/rt_  # w_e = k_ * v
# (a - b k v) k = Crr m g + ca v^2   ->  ca v^2 + b k^2 v + Crr m g - a k = 0
A_, B_, C_ = ca_, b_*k_**2, Crr_*m_*9.81 - a_*k_
v_ss = (-B_ + np.sqrt(B_**2 - 4*A_*C_))/(2*A_)
print(f"{'drive_ss_speed_mps':32s} {v_ss:.6g}")
print(f"{'drive_ss_engine_rpm':32s} {k_*v_ss*60/(2*np.pi):.6g}")
# --- shaft: simply supported L=0.3, d=0.025, P=1000 at midspan, steel
E, nu, d, Ls, P = 200e9, 0.3, 0.025, 0.3, 1000.0
I = np.pi*d**4/64; A = np.pi*d**2/4; G = E/(2*(1+nu))
kap = 6*(1+nu)/(7+6*nu)  # Cowper solid circle
print(f"{'shaft_ss_defl_EB':32s} {P*Ls**3/(48*E*I):.6g}")
print(f"{'shaft_ss_defl_timo':32s} {P*Ls**3/(48*E*I)+P*Ls/(4*kap*G*A):.6g}")
print(f"{'shaft_cowper_kappa':32s} {kap:.6g}")
# overhung: bearings at z=0, z=0.2, load P at z=0.25
a2, c2 = 0.2, 0.05
print(f"{'overhung_R1':32s} {-P*c2/a2:.6g}")
print(f"{'overhung_R2':32s} {P*(a2+c2)/a2:.6g}")
print(f"{'overhung_tip_defl_EB':32s} {P*c2**2*(a2+c2)/(3*E*I):.6g}")
# --- FE cantilever  L=10,h=w=1,E=1000,nu=0,P=1
print(f"{'fe_cant_tip':32s} {1*10**3/(3*1000/12) + 10/((5/6)*500*1):.6g}")
# string discrete N=20 segments, T0=100 N, m'=0.345 kg/m, span 0.2 m
N, T0, mp, sp_ = 20, 100.0, 0.345, 0.2
l_ = sp_/N; mnode = mp*l_/(1 + T0/EA)  # mass fixed by unstretched length
f_disc = (1/np.pi)*np.sqrt(T0/(mnode*l_))*np.sin(np.pi/(2*N))
print(f"{'string_f1_discrete':32s} {f_disc:.6g}")
print(f"{'string_f1_continuum':32s} {np.sqrt(T0/mp)/(2*sp_):.6g}")
# Kt shoulder fits (UNVERIFIED coeffs) at D/d=1.2 r/d=0.1 -> t=0.1d, r=0.1d
def kt_bend(D, d, r):
    t=(D-d)/2; q=t/r; s=np.sqrt(q)
    C=[0.947+1.206*s-0.131*q, 0.022-3.405*s+0.915*q, 0.869+1.777*s-0.555*q, -0.810+0.422*s-0.260*q]
    x=2*t/D; return C[0]+C[1]*x+C[2]*x**2+C[3]*x**3
print(f"{'kt_bend_D1.2_r0.1':32s} {kt_bend(1.2,1.0,0.1):.6g}")

# --- 09: profiled ramp and helix, vehicle with tire slip, belt heating
from numpy import arctanh
n9, mw9, rw9 = 3, 0.275 + 0.0387, 0.04445
w3k = 3000*2*np.pi/60
u9 = 0.01  # T9.2 curved ramp r_w = r_w0 + 0.4 u + 4 u^2, kappa = 0.4 + 8 u
print(f"{'t9_2_F_fw_u10mm_3000rpm':32s} {n9*mw9*(rw9 + 0.4*u9 + 4*u9**2)*(0.4 + 8*u9)*w3k**2:.17g}")
rc9 = 0.0381  # T9.4 helix: 45 deg to 10 mm, then 35 deg; h' = 1/(r_c tan alpha)
hp45, hp35 = 1/(rc9*np.tan(np.radians(45))), 1/(rc9*np.tan(np.radians(35)))
h9 = lambda x: x*hp45 if x <= 0.01 else 0.01*hp45 + (x - 0.01)*hp35
for xmm in (5, 15, 20):
    print(f"{f't9_4_h_{xmm}mm':32s} {h9(xmm/1e3):.17g}")
print(f"{'t9_4_hprime_15mm':32s} {hp35:.17g}")
Fl0_9 = 435.63 - 3.930*np.radians(30)/(rc9*np.tan(np.radians(45)))
print(f"{'t9_4_Fs_15mm_20Nm':32s} {Fl0_9 + 3353.66*0.015 + (20/2 + 3.930*(np.radians(30) + h9(0.015)))*hp35:.17g}")
# T9.6 fixed ratio, linear engine T = a - b w; dirt Crr, CdA 0.9, eta 0.9
a9, b9, R9, G9, rt9, m9, g9 = 40.0, 0.05, 1.0, 8.1, 0.2921, 250.0, 9.81
Crr9, ca9, eta9, mut9, s09, bL9 = 0.06, 0.5*1.2*0.9, 0.90, 0.7, 0.08, 0.60
k9 = R9*G9/rt9
A9, B9, C9 = ca9, eta9*b9*k9**2, Crr9*m9*g9 - eta9*a9*k9
print(f"{'t9_6_v_quadratic_noslip':32s} {(-B9 + np.sqrt(B9**2 - 4*A9*C9))/(2*A9):.17g}")
# with the tire: wheel speed v (1+s), s = s0 artanh(F/(mu_t N_r)), N_r = m g b/L at steady state
F9 = lambda v: Crr9*m9*g9 + ca9*v**2
s9 = lambda v: s09*arctanh(F9(v)/(mut9*m9*g9*bL9))
v9 = brentq(lambda v: eta9*k9*(a9 - b9*k9*v*(1 + s9(v))) - F9(v), 1.0, 30.0, xtol=1e-14)
print(f"{'t9_6_v_with_slip':32s} {v9:.17g}")
print(f"{'t9_6_tire_slip':32s} {s9(v9):.17g}")
# T9.7 tire on a dyno at v = 10 m/s, axle torque 200 N m: F_x = T/r_t, s = s0 artanh(F_x/(mu_t m g b/L))
print(f"{'t9_7_slip_T200_v10':32s} {s09*arctanh(200/rt9/(mut9*m9*g9*bL9)):.17g}")
print(f"{'t9_7_traction_limit_N':32s} {mut9*m9*g9*bL9:.17g}")
# T9.8 belt heating: m = rho A L, C_b = m c_p, A_s = L (b + 2 h_c), tau = C_b/(h A_s)
Lb9 = 0.895  # belt pitch length, dr26_default.yaml
Cb9 = 1150.0*3.0e-4*Lb9*1700.0
hA9 = 40.0*Lb9*(0.030 + 2*0.012)
print(f"{'t9_8_C_b':32s} {Cb9:.17g}")
print(f"{'t9_8_hA':32s} {hA9:.17g}")
print(f"{'t9_8_tau':32s} {Cb9/hA9:.17g}")
print(f"{'t9_8_rise_fraction_at_tau':32s} {1 - np.exp(-1):.17g}")

# --- 11: firing pulses, PTO torsion, clutch friction (11 sections 1-3)
from scipy.integrate import quad
B11, S11, rc11, l11, n11, p11 = 0.089, 0.069, 8.0, 0.105, 1.30, 101325.0
r11 = S11/2; A11 = np.pi*B11**2/4; Vc11 = A11*S11/(rc11 - 1); V111 = Vc11 + A11*S11
ths11, dth11, a11, m11 = np.radians(-15.0), np.radians(50.0), 5.0, 2.0
s11 = lambda th: r11*(1 - np.cos(th)) + l11*(1 - np.sqrt(1 - (r11/l11)**2*np.sin(th)**2))
print(f"{'t11_1_s_180':32s} {s11(np.pi):.17g}")
# T11.2: governed map mean at chart points (rpm_cut 3800, band 100 -> factor 1 up to 3700 rpm)
for rpm, T in ((2400, 25.08), (3000, 22.51), (3600, 18.30)):
    print(f"{f't11_2_Tmap_{rpm}':32s} {T*min(1.0, max(0.0, (3800 - rpm)/100)):.17g}")
# gas-torque cycle mean per joule of heat: first law over the closed part,
# W_Q = int x_b' (1 - (V/V1)^(n-1)) dpsi, mean torque = W_Q/(4 pi)
def dxb11(ps):
    t = max(ps - ths11, 0.0)/dth11
    return a11*(m11 + 1)*t**m11*np.exp(-a11*t**(m11 + 1))/dth11
WQ11 = quad(lambda ps: dxb11(ps)*(1 - ((Vc11 + A11*s11(ps))/V111)**(n11 - 1)), ths11, np.pi,
            points=[ths11 + k*dth11 for k in (0.5, 1, 1.5, 2) if ths11 + k*dth11 < np.pi],
            limit=400, epsabs=0.0, epsrel=1e-13)[0]
print(f"{'t11_gas_mean_per_J':32s} {WQ11/(4*np.pi):.17g}")
print(f"{'t11_heat_per_cycle_25Nm':32s} {4*np.pi*25.08/WQ11:.17g}")
# T11.5: PTO k = G pi d_eff^4/(32 L_eff), d_eff = d - w_key/2; two-inertia frequency
G11 = 200e9/(2*1.3); d11 = 0.0254 - 0.00635/2; L11 = 0.0886 + 0.030
k11 = G11*np.pi*d11**4/(32*L11)
Jf11, Jc11 = 0.035, 0.015
print(f"{'t11_5_k_pto':32s} {k11:.17g}")
print(f"{'t11_5_f_spec':32s} {np.sqrt(k11*(1/Jf11 + 1/Jc11))/(2*np.pi):.17g}")
# free-vibration test case: primary open at x_p = -g (8.9662 mm before contact), the three
# flyweights (0.3137 kg) at r_w = r_w0 - kappa g add to the clutch inertia
Jc11_t = Jc11 + 3*(0.275 + 0.0387)*(0.04445 - 0.5*0.0089662)**2
print(f"{'t11_5_f_open_clutch':32s} {np.sqrt(k11*(1/Jf11 + 1/Jc11_t))/(2*np.pi):.17g}")
# T11.6: 45 deg helix, tau_mov = 20 N m at x = 0: tau_c = tau_mov + k_t theta_t0
tauc11 = 20.0 + 3.930*np.radians(30.0); mu11 = 0.10
print(f"{'t11_6_F_close_opening':32s} {tauc11/0.0381*(1 + 2*mu11):.17g}")
print(f"{'t11_6_F_close_closing':32s} {tauc11/0.0381*(1 - 2*mu11):.17g}")
print(f"{'t11_6_ratio':32s} {(1 + 2*mu11)/(1 - 2*mu11):.17g}")
