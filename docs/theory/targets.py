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
# geometry from ratio range, strokes, beta
sp, ss = 25.48e-3, 24.48e-3
dp, ds = sp/(2*np.tan(beta)), ss/(2*np.tan(beta))
a = (ds + 0.8*dp)/(3-0.8)
rp = (a, a+dp); rs = (3*a, 0.8*(a+dp))
out.update(rp_min=rp[0], rp_max=rp[1], rs_max=rs[0], rs_min=rs[1])
def Lbelt(r1, r2, C):
    ph = np.arcsin((r2-r1)/C)
    return 2*C*np.cos(ph) + np.pi*(r1+r2) + 2*ph*(r2-r1)
out["note_Ldiff_vs_C_200mm"] = Lbelt(rp[0], rs[0], 0.2) - Lbelt(rp[1], rs[1], 0.2)
for k_, v in out.items():
    print(f"{k_:32s} {v:.6g}")

# --- derived center distance: equal belt length at both ratio extremes
Cc = brentq(lambda C: Lbelt(rp[0], rs[0], C) - Lbelt(rp[1], rs[1], C), 0.05, 2.0)
print(f"{'center_distance_C':32s} {Cc:.6g}")
print(f"{'belt_length_L':32s} {Lbelt(rp[0], rs[0], Cc):.6g}")
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
