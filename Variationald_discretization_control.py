import numpy as np
import math
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path

# ============================================================
# P0 constant control discretization vs variational control discretization
# Problem: -y'' + u y = f, y(0)=y(1)=0
# PLUS-sign projection: u = Proj_[a,b]((1/alpha) y p)
# ============================================================

alpha = 0.01
a, b = 0.1, 0.3
N_list = [4, 8, 16, 32, 64, 128, 256]
max_iter_constant = 300
max_iter_variational = 300
variational_abs_tol = 1e-10
variational_rel_tol = 1e-8
fixed_point_tol = 1e-12
out_dir = Path("merged_outputs")
out_dir.mkdir(exist_ok=True)

# ============================================================
# Exact data
# ============================================================
def y_exact(x):
    return x * (1.0 - x)

def p_exact(x):
    return 0.05 * np.sin(np.pi * x)

def proj(v):
    return np.minimum(b, np.maximum(a, v))

def u_exact(x):
    return proj((1.0 / alpha) * y_exact(x) * p_exact(x))

def f_fun(x):
    return 2.0 + u_exact(x) * y_exact(x)

def y_d_fun(x):
    minus_p_dd = 0.05 * (np.pi ** 2) * np.sin(np.pi * x)
    return y_exact(x) - (minus_p_dd + u_exact(x) * p_exact(x))

# ============================================================
# Quadrature and utilities
# ============================================================
GL_xi = np.array([-0.8611363116, -0.3399810436, 0.3399810436, 0.8611363116])
GL_w = np.array([0.3478548451, 0.6521451549, 0.6521451549, 0.3478548451])

def gauss_on_element(L, R):
    mid = 0.5 * (L + R)
    half = 0.5 * (R - L)
    return mid + half * GL_xi, half * GL_w

def compute_eoc(err, h, eps=1e-15):
    err = np.asarray(err, dtype=float)
    h = np.asarray(h, dtype=float)
    eoc = np.full_like(err, np.nan, dtype=float)
    for i in range(1, len(err)):
        if err[i] >= eps and err[i - 1] >= eps:
            eoc[i] = np.log(err[i] / err[i - 1]) / np.log(h[i] / h[i - 1])
    return eoc

def ref_line(hs, errs, order):
    return errs[0] * (hs / hs[0]) ** order

def l2_grid_error(v_exact_vals, v_h_vals, xg):
    return float(np.sqrt(np.trapezoid((v_exact_vals - v_h_vals) ** 2, xg)))

def eval_p1_from_interior(N, V_int, xg):
    x_nodes = np.linspace(0.0, 1.0, N + 1)
    V = np.zeros(N + 1)
    V[1:N] = V_int
    return eval_p1_from_nodes(x_nodes, V, xg)

def eval_p1_from_nodes(x_nodes, V_nodes, xg):
    xg = np.asarray(xg)
    vg = np.zeros_like(xg, dtype=float)
    N = len(x_nodes) - 1
    h = x_nodes[1] - x_nodes[0]
    mask = (xg > x_nodes[0]) & (xg < x_nodes[-1])
    xm = xg[mask]
    e = np.floor(xm / h).astype(int)
    e = np.clip(e, 0, N - 1)
    L = x_nodes[e]
    R = x_nodes[e + 1]
    VL = V_nodes[e]
    VR = V_nodes[e + 1]
    vg[mask] = VL + (VR - VL) * (xm - L) / (R - L)
    return vg

# ============================================================
# Constant/P0 control solver
# ============================================================
def assemble_K_M0_constant(N):
    h = 1.0 / N
    n = N - 1
    mainK = (2.0 / h) * np.ones(n)
    offK = (-1.0 / h) * np.ones(n - 1)
    K = np.diag(mainK) + np.diag(offK, 1) + np.diag(offK, -1)
    mainM = (4.0 * h / 6.0) * np.ones(n)
    offM = (h / 6.0) * np.ones(n - 1)
    M0 = np.diag(mainM) + np.diag(offM, 1) + np.diag(offM, -1)
    return K, M0

def assemble_reaction_Mu_constant(N, u_elem):
    h = 1.0 / N
    n = N - 1
    Mu = np.zeros((n, n))
    for e in range(1, N + 1):
        ue = u_elem[e - 1]
        Me = ue * (h / 6.0) * np.array([[2.0, 1.0], [1.0, 2.0]])
        gL, gR = e - 1, e
        for li, gi in enumerate([gL, gR]):
            for lj, gj in enumerate([gL, gR]):
                if 1 <= gi <= N - 1 and 1 <= gj <= N - 1:
                    Mu[gi - 1, gj - 1] += Me[li, lj]
    return Mu

def assemble_load_vector_constant(N, gfun):
    x_nodes = np.linspace(0.0, 1.0, N + 1)
    n = N - 1
    Lvec = np.zeros(n)
    for e in range(1, N + 1):
        L, R = x_nodes[e - 1], x_nodes[e]
        h = R - L
        xq, wq = gauss_on_element(L, R)
        phiL = (R - xq) / h
        phiR = (xq - L) / h
        gq = gfun(xq)
        gL, gR = e - 1, e
        if 1 <= gL <= N - 1:
            Lvec[gL - 1] += np.sum(wq * gq * phiL)
        if 1 <= gR <= N - 1:
            Lvec[gR - 1] += np.sum(wq * gq * phiR)
    return Lvec

def update_control_constant(N, Y_int, P_int):
    h = 1.0 / N
    Y = np.zeros(N + 1)
    P = np.zeros(N + 1)
    Y[1:N] = Y_int
    P[1:N] = P_int
    u_new = np.zeros(N)
    for e in range(1, N + 1):
        YL, YR = Y[e - 1], Y[e]
        PL, PR = P[e - 1], P[e]
        integral_yp = (h / 6.0) * (2 * YL * PL + YL * PR + YR * PL + 2 * YR * PR)
        u_new[e - 1] = float(proj((1.0 / (alpha * h)) * integral_yp))
    return u_new

def solve_constant_N(N, u0=0.2):
    K, M0 = assemble_K_M0_constant(N)
    F = assemble_load_vector_constant(N, f_fun)
    Yd = assemble_load_vector_constant(N, y_d_fun)
    u = np.full(N, u0, dtype=float)
    it_count = max_iter_constant
    for it in range(max_iter_constant):
        A = K + assemble_reaction_Mu_constant(N, u)
        Y = np.linalg.solve(A, F)
        P = np.linalg.solve(A, M0 @ Y - Yd)
        u_new = update_control_constant(N, Y, P)
        if np.max(np.abs(u_new - u)) < fixed_point_tol:
            u = u_new
            it_count = it + 1
            break
        u = u_new
    A = K + assemble_reaction_Mu_constant(N, u)
    Y = np.linalg.solve(A, F)
    P = np.linalg.solve(A, M0 @ Y - Yd)
    return u, Y, P, it_count

def eval_constant_control(N, u_elem, xg):
    xg = np.asarray(xg)
    h = 1.0 / N
    idx = np.floor(xg / h).astype(int)
    idx = np.clip(idx, 0, N - 1)
    return u_elem[idx]

def l2_control_error_constant(N, u_elem):
    x_nodes = np.linspace(0.0, 1.0, N + 1)
    err2 = 0.0
    for e in range(1, N + 1):
        L, R = x_nodes[e - 1], x_nodes[e]
        xq, wq = gauss_on_element(L, R)
        diff = u_elem[e - 1] - u_exact(xq)
        err2 += np.sum(wq * diff * diff)
    return float(np.sqrt(err2))

# ============================================================
# Variational/segmented control solver
# ============================================================
def poly_mult(p, q):
    res = [0.0] * (len(p) + len(q) - 1)
    for i, pi in enumerate(p):
        for j, qj in enumerate(q):
            res[i + j] += pi * qj
    return res

def poly_integral_on_segment(coeffs, xs, xe):
    return sum(ck * (xe ** (k + 1) - xs ** (k + 1)) / (k + 1) for k, ck in enumerate(coeffs))

def project_control_on_element(xL, xR, A, B, C):
    def q(x):
        return A * x * x + B * x + C
    values = [q(xL), q(xR)]
    if abs(A) > 1e-14:
        xV = -B / (2 * A)
        if xL <= xV <= xR:
            values.append(q(xV))
    qmin, qmax = min(values), max(values)
    if qmax <= a:
        return [(xL, xR, 'const_a', None)]
    if qmin >= b:
        return [(xL, xR, 'const_b', None)]
    if qmin >= a and qmax <= b:
        return [(xL, xR, 'quad', (A, B, C))]
    roots = []
    def add_roots(level):
        AA, BB, CC = A, B, C - level
        if abs(AA) < 1e-14:
            if abs(BB) > 1e-14:
                r = -CC / BB
                if xL <= r <= xR:
                    roots.append(r)
            return
        disc = BB * BB - 4 * AA * CC
        if disc < -1e-14:
            return
        disc = max(0.0, disc)
        sd = math.sqrt(disc)
        for r in [(-BB - sd) / (2 * AA), (-BB + sd) / (2 * AA)]:
            if xL <= r <= xR:
                roots.append(r)
    if qmin < a < qmax:
        add_roots(a)
    if qmin < b < qmax:
        add_roots(b)
    pts = sorted(set([xL, xR] + roots))
    segments = []
    for i in range(len(pts) - 1):
        xs, xe = pts[i], pts[i + 1]
        qmid = q(0.5 * (xs + xe))
        if qmid <= a:
            segments.append((xs, xe, 'const_a', None))
        elif qmid >= b:
            segments.append((xs, xe, 'const_b', None))
        else:
            segments.append((xs, xe, 'quad', (A, B, C)))
    return segments

def local_C_from_segments(xL, xR, segments):
    h = xR - xL
    phi1 = [xR / h, -1.0 / h]
    phi2 = [-xL / h, 1.0 / h]
    C_loc = np.zeros((2, 2))
    for xs, xe, typ, data in segments:
        if typ == 'const_a':
            u_poly = [a]
        elif typ == 'const_b':
            u_poly = [b]
        else:
            A, B, C = data
            u_poly = [C, B, A]
        for m, pm in enumerate([phi1, phi2]):
            for n, pn in enumerate([phi1, phi2]):
                C_loc[m, n] += poly_integral_on_segment(poly_mult(u_poly, poly_mult(pm, pn)), xs, xe)
    return C_loc

def assemble_C_global_variational(x, all_segments):
    N = len(x) - 1
    C = np.zeros((N + 1, N + 1))
    for e in range(N):
        C_loc = local_C_from_segments(x[e], x[e + 1], all_segments[e])
        C[e:e + 2, e:e + 2] += C_loc
    return C

def qp_coefficients_on_element(xL, xR, Yi, Yip1, Pi, Pip1):
    h = xR - xL
    sy = (Yip1 - Yi) / h
    cy = Yi - sy * xL
    sp = (Pip1 - Pi) / h
    cp = Pi - sp * xL
    return (sy * sp / alpha, (sy * cp + cy * sp) / alpha, cy * cp / alpha)

def assemble_K_M_F_variational(x, f_handle):
    N = len(x) - 1
    K = np.zeros((N + 1, N + 1))
    M = np.zeros((N + 1, N + 1))
    F = np.zeros(N + 1)
    for e in range(N):
        xL, xR = x[e], x[e + 1]
        h = xR - xL
        K_loc = (1.0 / h) * np.array([[1, -1], [-1, 1]])
        M_loc = (h / 6.0) * np.array([[2, 1], [1, 2]])
        xq, wq = gauss_on_element(xL, xR)
        phi1 = (xR - xq) / h
        phi2 = (xq - xL) / h
        fq = f_handle(xq)
        F_loc = np.array([np.sum(wq * fq * phi1), np.sum(wq * fq * phi2)])
        K[e:e + 2, e:e + 2] += K_loc
        M[e:e + 2, e:e + 2] += M_loc
        F[e:e + 2] += F_loc
    return K, M, F

def solve_variational_N(N):
    x = np.linspace(0.0, 1.0, N + 1)
    K, M, F = assemble_K_M_F_variational(x, f_fun)

    Yd_nodes = y_d_fun(x)

    I = np.arange(1, N)
    K_int = K[np.ix_(I, I)]
    M_int = M[np.ix_(I, I)]
    F_int = F[I]
    Yd_int = Yd_nodes[I]

    segments = [[(x[e], x[e + 1], 'const_b', None)] for e in range(N)]

    Y = np.zeros(N + 1)
    P = np.zeros(N + 1)

    it_count = max_iter_variational
    converged = False
    control_change_L2 = np.nan
    relative_change = np.nan

    for it in range(max_iter_variational):
        C = assemble_C_global_variational(x, segments)
        A_int = K_int + C[np.ix_(I, I)]

        Y_int = np.linalg.solve(A_int, F_int)
        Y[:] = 0.0
        Y[I] = Y_int

        P_int = np.linalg.solve(A_int, M_int @ (Y_int - Yd_int))
        P[:] = 0.0
        P[I] = P_int

        new_segments = []

        for e in range(N):
            Aq, Bq, Cq = qp_coefficients_on_element(
                x[e], x[e + 1],
                Y[e], Y[e + 1],
                P[e], P[e + 1]
            )

            new_segments.append(
                project_control_on_element(x[e], x[e + 1], Aq, Bq, Cq)
            )

        # ====================================================
        # tolerance-based stopping rule for variational control
        # ====================================================
        Xtest = np.linspace(0.0, 1.0, max(4000, 20 * N))

        u_old = eval_variational_control(x, segments, Xtest)
        u_new = eval_variational_control(x, new_segments, Xtest)

        control_change_L2 = float(
            np.sqrt(np.trapezoid((u_new - u_old) ** 2, Xtest))
        )

        norm_new = float(
            np.sqrt(np.trapezoid(u_new ** 2, Xtest))
        )

        if control_change_L2 <= variational_abs_tol + variational_rel_tol * max(1.0, norm_new):
            segments = new_segments
            it_count = it + 1
            #converged = True
            break

        segments = new_segments

    return x, segments, Y, P, K, M, it_count, control_change_L2

def eval_variational_control(x, segments, Xfine):
    Q = np.zeros_like(Xfine, dtype=float)
    for e in range(len(segments)):
        xL, xR = x[e], x[e + 1]
        if e == len(segments) - 1:
            mask = (Xfine >= xL) & (Xfine <= xR)
        else:
            mask = (Xfine >= xL) & (Xfine < xR)
        Xloc = Xfine[mask]
        vals = np.zeros_like(Xloc)
        for xs, xe, typ, data in segments[e]:
            seg_mask = (Xloc >= xs - 1e-14) & (Xloc <= xe + 1e-14)
            if typ == 'const_a':
                vals[seg_mask] = a
            elif typ == 'const_b':
                vals[seg_mask] = b
            else:
                A, B, C = data
                vals[seg_mask] = A * Xloc[seg_mask] ** 2 + B * Xloc[seg_mask] + C
        Q[mask] = vals
    return Q

def l2_control_error_variational(x, segments):
    err = 0.0
    for e in range(len(segments)):
        xL, xR = x[e], x[e + 1]
        xq, wq = gauss_on_element(xL, xR)
        uh = eval_variational_control(x, segments, xq)
        diff = uh - u_exact(xq)
        err += np.sum(wq * diff * diff)
    return float(np.sqrt(err))

# ============================================================
# Run studies
# ============================================================
def run_constant_study(x_plot):
    rows = []
    solutions = {}
    y_ex_vals = y_exact(x_plot)
    p_ex_vals = p_exact(x_plot)
    for N in N_list:
        uh, Y, P, iters = solve_constant_N(N)
        yh = eval_p1_from_interior(N, Y, x_plot)
        ph = eval_p1_from_interior(N, P, x_plot)
        rows.append({
            'method': 'constant_P0', 'N': N, 'h': 1.0 / N, 'iters': iters,
            'err_y_L2': l2_grid_error(y_ex_vals, yh, x_plot),
            'err_p_L2': l2_grid_error(p_ex_vals, ph, x_plot),
            'err_u_L2': l2_control_error_constant(N, uh),
        })
        solutions[N] = {'u': uh, 'Y_int': Y, 'P_int': P}
    df = pd.DataFrame(rows)
    for name in ['y', 'p', 'u']:
        df[f'EOC_{name}'] = compute_eoc(df[f'err_{name}_L2'].values, df['h'].values)
    return df, solutions

def run_variational_study(x_plot):
    rows = []
    solutions = {}
    y_ex_vals = y_exact(x_plot)
    p_ex_vals = p_exact(x_plot)
    for N in N_list:
        x, segments, Y, P, K, M, iters, control_change_L2 = solve_variational_N(N)
        yh = eval_p1_from_nodes(x, Y, x_plot)
        ph = eval_p1_from_nodes(x, P, x_plot)
        rows.append({
            'method': 'variational_segmented', 'N': N, 'h': 1.0 / N, 'iters': iters,
            'control_change_L2': control_change_L2,
            'err_y_L2': l2_grid_error(y_ex_vals, yh, x_plot),
            'err_p_L2': l2_grid_error(p_ex_vals, ph, x_plot),
            'err_u_L2': l2_control_error_variational(x, segments),
        })
        solutions[N] = {'x': x, 'segments': segments, 'Y': Y, 'P': P}
    df = pd.DataFrame(rows)
    for name in ['y', 'p', 'u']:
        df[f'EOC_{name}'] = compute_eoc(df[f'err_{name}_L2'].values, df['h'].values)
    return df, solutions

def save_tables(df_constant, df_variational):
    combined = pd.concat([df_constant, df_variational], ignore_index=True)
    df_constant.to_csv(out_dir / 'eoc_constant_control.csv', index=False)
    df_variational.to_csv(out_dir / 'eoc_variational_control.csv', index=False)
    combined.to_csv(out_dir / 'eoc_combined_constant_vs_variational.csv', index=False)
    with pd.ExcelWriter(out_dir / 'eoc_tables_constant_vs_variational.xlsx') as writer:
        df_constant.to_excel(writer, sheet_name='constant_P0', index=False)
        df_variational.to_excel(writer, sheet_name='variational', index=False)
        combined.to_excel(writer, sheet_name='combined', index=False)
    return combined

def make_plots(df_c, df_v, sol_c, sol_v, x_plot):
    y_ex_vals = y_exact(x_plot)
    p_ex_vals = p_exact(x_plot)
    u_ex_vals = u_exact(x_plot)

    # 1) EOC/error comparison: constant vs variational for control
    plt.figure(figsize=(8, 5))
    plt.loglog(df_c['h'], df_c['err_u_L2'], marker='o', label='constant/P0 control error')
    plt.loglog(df_v['h'], df_v['err_u_L2'], marker='s', label='variational segmented control error')
    plt.loglog(df_c['h'], ref_line(df_c['h'].values, df_c['err_u_L2'].values, 1), '--', label='ref O(h)')
    plt.loglog(df_v['h'], ref_line(df_v['h'].values, df_v['err_u_L2'].values, 2), '--', label='ref O(h²)')
    plt.gca().invert_xaxis()
    plt.xlabel('h = 1/N')
    plt.ylabel(r'$L^2$ control error')
    plt.title('Control convergence: constant vs variational')
    plt.grid(True, which='both')
    plt.legend()
    plt.tight_layout()
    #plt.savefig(out_dir / 'eoc_control_constant_vs_variational.png', dpi=200)
    plt.show()

    # 2) All variables comparison by method
    # for variable, err_col, title in [('state', 'err_y_L2', 'State'), ('adjoint', 'err_p_L2', 'Adjoint'), ('control', 'err_u_L2', 'Control')]:
    #     plt.figure(figsize=(8, 5))
    #     plt.loglog(df_c['h'], df_c[err_col], marker='o', label=f'constant/P0 {variable}')
    #     plt.loglog(df_v['h'], df_v[err_col], marker='s', label=f'variational {variable}')
    #     plt.gca().invert_xaxis()
    #     plt.xlabel('h = 1/N')
    #     plt.ylabel(r'$L^2$ error')
    #     plt.title(f'{title} EOC comparison: constant vs variational')
    #     plt.grid(True, which='both')
    #     plt.legend()
    #     plt.tight_layout()
    #     plt.savefig(out_dir / f'eoc_{variable}_constant_vs_variational.png', dpi=200)
    #     plt.show()

    # 3) Exact vs both numerical methods on finest mesh
    Nf = N_list[-1]
    c = sol_c[Nf]
    v = sol_v[Nf]
    y_c = eval_p1_from_interior(Nf, c['Y_int'], x_plot)
    p_c = eval_p1_from_interior(Nf, c['P_int'], x_plot)
    u_c = eval_constant_control(Nf, c['u'], x_plot)
    y_v = eval_p1_from_nodes(v['x'], v['Y'], x_plot)
    p_v = eval_p1_from_nodes(v['x'], v['P'], x_plot)
    u_v = eval_variational_control(v['x'], v['segments'], x_plot)

    plot_data = [
        ('state', y_ex_vals, y_v, y_c, 'y'),
        ('adjoint', p_ex_vals, p_v, p_c, 'p'),
        ('control', u_ex_vals, u_v, u_c, 'u'),
    ]
    for name, exact_vals, var_vals, const_vals, ylabel in plot_data:
        plt.figure(figsize=(10, 5))
        plt.plot(x_plot, exact_vals, label=f'exact {ylabel}', linewidth=2)
        plt.plot(x_plot, var_vals, '--', label=f'variational {ylabel}_h, N={Nf}')
        plt.plot(x_plot, const_vals, ':', label=f'constant/P0 {ylabel}_h, N={Nf}')
        if name == 'control':
            plt.hlines([a, b], 0.0, 1.0, linestyles='dashed', label='bounds')
        plt.xlabel('x')
        plt.ylabel(ylabel)
        plt.title(f'{name.capitalize()}: exact vs variational vs constant')
        plt.grid(True)
        plt.legend()
        plt.tight_layout()
        #plt.savefig(out_dir / f'exact_vs_variational_vs_constant_{name}.png', dpi=200)
        plt.show()


def print_table(df, title):
    display = df.copy()
    display['h'] = display['h'].map(lambda z: f'{z:.6f}')
    for col in ['err_y_L2', 'err_p_L2', 'err_u_L2']:
        display[col] = display[col].map(lambda z: f'{z:.6e}')
    for col in ['EOC_y', 'EOC_p', 'EOC_u']:
        display[col] = display[col].map(lambda z: '-' if pd.isna(z) else f'{z:.4f}')
    print(f'\n{title}')
    print(display.to_string(index=False))


def main():
    x_plot = np.linspace(0.0, 1.0, 4001)
    df_constant, sol_constant = run_constant_study(x_plot)
    df_variational, sol_variational = run_variational_study(x_plot)
    #combined = save_tables(df_constant, df_variational)
    print_table(df_constant, 'EOC table: constant/P0 control')
    print_table(df_variational, 'EOC table: variational segmented control')
    make_plots(df_constant, df_variational, sol_constant, sol_variational, x_plot)
    # print(f'\nSaved outputs in: {out_dir.resolve()}')
    # print('Files:')
    # for path in sorted(out_dir.iterdir()):
    #     print(f'  - {path.name}')

if __name__ == '__main__':
    main()