import pandas as pd
import numpy as np
import datetime
import mpmath
from sympy.physics.wigner import wigner_3j
from tqdm import tqdm

# Set arbitrary precision to 50 decimal places to prevent catastrophic cancellation
mpmath.mp.dps = 50
mpmath.mp.maxprec = 50000
mpmath.mp.zeroprec = 200

def spherical_jn_mp(k, a):
    """Evaluates the spherical Bessel function of the first kind using mpmath."""
    if a == 0:
        return mpmath.mpf(1) if k == 0 else mpmath.mpf(0)
    return mpmath.sqrt(mpmath.pi / (2 * a)) * mpmath.besselj(k + 0.5, a)

def calculate_A1_integral(l, lp, m, mp, a, k_max=None, tol=1e-11, max_k_safety=250):
    scale_factor = mpmath.mpf(1.0)
    if m < 0: scale_factor *= ((-1)**abs(m)) * mpmath.fac(l - abs(m)) / mpmath.fac(l + abs(m))
    if mp < 0: scale_factor *= ((-1)**abs(mp)) * mpmath.fac(lp - abs(mp)) / mpmath.fac(lp + abs(mp))
        
    m_abs, mp_abs = abs(m), abs(mp)
    M = m_abs + mp_abs
    Lambda = mpmath.sqrt((mpmath.fac(lp + mp_abs) * mpmath.fac(l + m_abs)) / (mpmath.fac(lp - mp_abs) * mpmath.fac(l - m_abs)))
    
    def C(k, j): return ((-1)**j * mpmath.fac(2*k - 2*j)) / ((2**k) * mpmath.fac(j) * mpmath.fac(k - j) * mpmath.fac(k - 2*j))
    def D(n, M, s): return ((-1)**s * mpmath.fac(2*n - 2*s)) / ((2**n) * mpmath.fac(s) * mpmath.fac(n - s) * mpmath.fac(n - M - 2*s))
    def B(p, M): return mpmath.mpf(0) if (p % 2 != 0 or p < 0) else (mpmath.gamma((p + 1) / 2.0) * mpmath.gamma((M + 2) / 2.0)) / mpmath.gamma((p + M + 3) / 2.0)

    # Determine turning-point cutoff horizon
    a_f = float(a)
    k_min_horizon = int(mpmath.ceil(a_f + 5.0 * (a_f**(1/3)) + 15.0)) if a_f > 0 else 0
    upper_bound = (k_max + 1) if k_max is not None else max_k_safety

    total_sum = mpmath.mpc(0, 0)
    consecutive_converged = 0

    for k in range(upper_bound):
        bessel_term = (mpmath.j**k) * (2*k + 1) * spherical_jn_mp(k, a)
        current_k_term = mpmath.mpc(0, 0)
        for n in range(abs(l - lp), l + lp + 1):
            w3j_1 = wigner_3j(l, lp, n, 0, 0, 0)
            w3j_2 = wigner_3j(l, lp, n, m_abs, mp_abs, -M)
            if w3j_1 == 0 or w3j_2 == 0: continue
            
            w1_mp = mpmath.mpf(str(w3j_1.evalf(50)))
            w2_mp = mpmath.mpf(str(w3j_2.evalf(50)))
            
            n_term = (2*n + 1) * w1_mp * w2_mp * mpmath.sqrt(mpmath.fac(n - M) / mpmath.fac(n + M))
            inner_sum = mpmath.mpf(0)
            for j in range(k // 2 + 1):
                c_val = C(k, j)
                for s in range((n - M) // 2 + 1):
                    p = k + n - M - 2*j - 2*s
                    inner_sum += c_val * D(n, M, s) * B(p, M)
            current_k_term += bessel_term * n_term * inner_sum

        total_sum += current_k_term

        # Check convergence once beyond turning-point horizon (when k_max is unspecified or dynamic)
        if k_max is None:
            term_magnitude = abs(scale_factor * Lambda * current_k_term)
            if k >= k_min_horizon and term_magnitude < tol:
                consecutive_converged += 1
                if consecutive_converged >= 3:
                    break
            else:
                consecutive_converged = 0

    return scale_factor * Lambda * total_sum

def calculate_A2_integral(l, lp, m, mp, a, k_max):
    scale_factor = mpmath.mpf(1.0)
    if m < 0: scale_factor *= ((-1)**abs(m)) * mpmath.fac(l - abs(m)) / mpmath.fac(l + abs(m))
    if mp < 0: scale_factor *= ((-1)**abs(mp)) * mpmath.fac(lp - abs(mp)) / mpmath.fac(lp + abs(mp))
        
    m_abs, mp_abs = abs(m), abs(mp)
    M = m_abs + mp_abs
    Lambda = mpmath.sqrt((mpmath.fac(lp + mp_abs) * mpmath.fac(l + m_abs)) / (mpmath.fac(lp - mp_abs) * mpmath.fac(l - m_abs)))

    def D(n, M, s): return ((-1)**s * mpmath.fac(2*n - 2*s)) / ((2**n) * mpmath.fac(s) * mpmath.fac(n - s) * mpmath.fac(n - M - 2*s))
    def B(p, M): return mpmath.mpf(0) if (p % 2 != 0 or p < 0) else (mpmath.gamma((p + 1) / 2.0) * mpmath.gamma((M + 2) / 2.0)) / mpmath.gamma((p + M + 3) / 2.0)

    def int_P_r_M(r, M):
        if (r - M) % 2 != 0 or r < M: return mpmath.mpf(0)
        return sum(D(r, M, s) * B(r - M - 2*s, M) for s in range((r - M) // 2 + 1))

    def C_M_func(n, M): return mpmath.sqrt(mpmath.fac(n + M) / mpmath.fac(n - M))
    def D_M_func(M, r): return ((-1)**M) * mpmath.sqrt(mpmath.fac(r - M) / mpmath.fac(r + M)) * int_P_r_M(r, M)

    upper_bound = (k_max + 1) if k_max is not None else 100
    total_sum = mpmath.mpc(0, 0)
    for k in range(upper_bound):
        bessel_term = (mpmath.j**k) * (2*k + 1) * spherical_jn_mp(k, a)
        for n in range(abs(l - lp), l + lp + 1):
            w3j_1 = wigner_3j(l, lp, n, 0, 0, 0)
            w3j_2 = wigner_3j(l, lp, n, m_abs, mp_abs, -M)
            if w3j_1 == 0 or w3j_2 == 0: continue
            
            w1_mp = mpmath.mpf(str(w3j_1.evalf(50)))
            w2_mp = mpmath.mpf(str(w3j_2.evalf(50)))
            
            n_term = (2*n + 1) * w1_mp * w2_mp * mpmath.sqrt(mpmath.fac(n - M) / mpmath.fac(n + M))
            c_m_val = C_M_func(n, M)
            
            inner_sum = mpmath.mpf(0)
            for r in range(abs(n - k), n + k + 1):
                w3j_r1 = wigner_3j(n, k, r, 0, 0, 0)
                w3j_r2 = wigner_3j(n, k, r, -M, 0, M) 
                if w3j_r1 == 0 or w3j_r2 == 0: continue
                
                wr1_mp = mpmath.mpf(str(w3j_r1.evalf(50)))
                wr2_mp = mpmath.mpf(str(w3j_r2.evalf(50)))
                    
                r_term = D_M_func(M, r) * (2*r + 1) * wr1_mp * wr2_mp
                inner_sum += r_term
                
            total_sum += bessel_term * n_term * c_m_val * inner_sum
    return scale_factor * Lambda * total_sum

def calculate_M0_exact_integral(l, lp, a):
    total_sum = mpmath.mpc(0, 0)
    for k in range(abs(l - lp), l + lp + 1):
        if (l + lp + k) % 2 != 0: continue
        w3j = wigner_3j(l, lp, k, 0, 0, 0)
        if w3j == 0: continue
        w_mp = mpmath.mpf(str(w3j.evalf(50)))
        bessel_term = (mpmath.j**k) * (2*k + 1) * spherical_jn_mp(k, a)
        total_sum += mpmath.mpf(2.0) * bessel_term * (w_mp**2)
    return total_sum

def numerical_ground_truth(l, lp, m, mp, a):
    # Parity check: P_l^m(-x) * P_lp^mp(-x) = (-1)^(l + lp + m + mp) * [P_l^m(x) * P_lp^mp(x)]
    prod_parity = (-1)**(l + lp + abs(m) + abs(mp))

    # Real part contains cos(ax) (even function).
    if prod_parity == -1:
        real_val = mpmath.mpf(0)
    else:
        def integrand_real(x):
            return mpmath.legenp(l, m, x) * mpmath.legenp(lp, mp, x) * mpmath.cos(a * x)
        real_val = 2 * mpmath.quad(integrand_real, [0, 1], method='tanh-sinh')

    # Imaginary part contains sin(ax) (odd function).
    if prod_parity == 1:
        imag_val = mpmath.mpf(0)
    else:
        def integrand_imag(x):
            return mpmath.legenp(l, m, x) * mpmath.legenp(lp, mp, x) * mpmath.sin(a * x)
        imag_val = 2 * mpmath.quad(integrand_imag, [0, 1], method='tanh-sinh')

    return mpmath.mpc(real_val, imag_val)

def explore_integrals(mode, l_in, lp_in, a_in, k_max=None, m_in=None, mp_in=None):
    """
    mode='full': Iterates all valid m, mp for given l, lp lists and a combinations.
    mode='scan': Sweeps across a_in array for a single state defined by integers (l, lp, m, mp, k_max).
                 Pass k_max=None to dynamically auto-converge to 1e-10 at each point.
    """
    results = []
    combinations = []

    a_vals = a_in if isinstance(a_in, (list, np.ndarray)) else [a_in]

    if mode == 'full':
        l_vals = l_in if isinstance(l_in, list) else [l_in]
        lp_vals = lp_in if isinstance(lp_in, list) else [lp_in]

        for l in l_vals:
            for lp in lp_vals:
                for a in a_vals:
                    for m in range(-l, l + 1):
                        for mp in range(-lp, lp + 1):
                            combinations.append((l, lp, m, mp, a))
        desc_text = "Evaluating Full State Grid"

    elif mode == 'scan':
        if any(v is None for v in (m_in, mp_in)):
            raise ValueError("Scan mode requires explicit m_in and mp_in integers.")

        l, lp, m, mp = l_in, lp_in, m_in, mp_in

        if abs(m) > l or abs(mp) > lp:
            print(f"Invalid state: m={m} or mp={mp} exceeds angular bounds.")
            return pd.DataFrame()

        for a in a_vals:
            combinations.append((l, lp, m, mp, a))

        desc_text = f"Scanning 'a' for |l={l}, m={m}> |lp={lp}, mp={mp}>"
    else:
        raise ValueError("Mode must be either 'full' or 'scan'.")

    for l, lp, m, mp, a in tqdm(combinations, desc=desc_text, unit="step"):
        try:
            a1_val = calculate_A1_integral(l, lp, m, mp, a, k_max=k_max)
            a2_val = calculate_A2_integral(l, lp, m, mp, a, k_max=k_max)
            num_val = numerical_ground_truth(l, lp, m, mp, a)

            # Convert to float and snap near-zero values
            vals = [
                float(a1_val.real), float(a1_val.imag), 
                float(a2_val.real), float(a2_val.imag), 
                float(num_val.real), float(num_val.imag)
            ]
            vals = [0.0 if abs(v) < 1e-12 else v for v in vals]
            a1_real, a1_imag, a2_real, a2_imag, num_real, num_imag = vals

            # Robust handling of M0 without assuming attributes exist when m != 0
            if m == 0 and mp == 0:
                m0_val = calculate_M0_exact_integral(l, lp, a)
                m0_real = 0.0 if abs(m0_val.real) < 1e-12 else float(m0_val.real)
                m0_imag = 0.0 if abs(m0_val.imag) < 1e-12 else float(m0_val.imag)
            else:
                m0_real, m0_imag = None, None

            results.append({
                'l': l, 'lp': lp, 'm': m, 'mp': mp, 'a': a,
                'A1 (Re)': a1_real, 'A2 (Re)': a2_real, 'M0 (Re)': m0_real, 'Num (Re)': num_real,
                'A1 (Im)': a1_imag, 'A2 (Im)': a2_imag, 'M0 (Im)': m0_imag, 'Num (Im)': num_imag,
            })
        except Exception as e:
            tqdm.write(f"Error evaluating l={l}, lp={lp}, m={m}, mp={mp}, a={a}: {e}")

    return pd.DataFrame(results)

if __name__ == "__main__":
    # --- CONFIGURATION ---
    active_mode = 'full'  # Toggle between 'full' or 'scan'
    date_str = datetime.datetime.now().strftime("%Y%m%d")

    if active_mode == 'full':
        l_in, lp_in = [0, 2], [0, 2]
        a_in = 13.2
        test_k_max = None

        df = explore_integrals('full', l_in, lp_in, a_in, k_max=test_k_max)

        def fmt(val): return f"{min(val)}-{max(val)}" if isinstance(val, list) else str(val)
        filename = f"full_l{fmt(l_in)}_lp{fmt(lp_in)}_a{fmt(a_in)}_kmax{test_k_max}_{date_str}.csv"

    elif active_mode == 'scan':
        # The 5 strict integer inputs for a targeted state
        scan_l = 1 
        scan_lp = 2 
        scan_m = 1
        scan_mp = 0 
        
        # Set to None for automatic a-dependent Wiscombe-Debye convergence to 1e-10,
        # or set to an integer (e.g. 1, 50, 100) to force a fixed cutoff.
        test_k_max = None 

        # The array being scanned
        a_array = np.linspace(0, 50, 251)

        df = explore_integrals('scan', scan_l, scan_lp, a_array, k_max=test_k_max, m_in=scan_m, mp_in=scan_mp)

        k_tag = f"kmax{test_k_max}" if test_k_max is not None else "kmaxAuto"
        filename = f"scan_l{scan_l}_m{scan_m}_lp{scan_lp}_mp{scan_mp}_{k_tag}_{date_str}.csv"

    # --- OUTPUT ---
    pd.set_option('display.max_columns', None)
    pd.set_option('display.width', 200)
    pd.set_option('display.float_format', '{:.10f}'.format)
    print(df.to_string(index=False, na_rep='-'))

    df.to_csv(filename, index=False, na_rep='-', float_format='%.10f')
    print(f"\nData exported to: {filename}")
