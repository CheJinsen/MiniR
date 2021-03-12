include("dpq.jl")
include("ref.jl")

function bratio(a::Real, b::Real, x::Real, y::Real, w::RefFloat, w1::RefFloat, ierr::RefInt, lp::Bool)
    do_swap = false
    n = 0
    z, a0, b0, x0, y0, λ = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    eps = 2.220446049250313e-16

    ierr1 = RefInt(0)

    w.val = D0(lp)
    w1.val = D0(lp)

    if isnan(x) || isnan(y) || isnan(a) || isnan(b)
        ierr.val = 9
        return
    end

    if a < 0.0 || b < 0.0
        ierr.val = 1
        return
    end
    if a == 0.0 && b == 0.0
        ierr.val = 2
        return
    end
    if x < 0.0 || x > 1.0
        ierr.val = 3
        return
    end
    if y < 0.0 || y > 1.0
        ierr.val = 4
        return
    end

    z = x + y - 1.0
    ierr.val = 0

    if abs(z) > eps * 3.0
        ierr.val = 5
        return
    end

    x == 0.0 && @goto L200
    y == 0.0 && @goto L210
    a == 0.0 && @goto L211
    b == 0.0 && @goto L201

    eps = max(eps, 1e-15)
    a_lt_b = a < b
    
    if (a_lt_b ? b : a) < eps * 0.001
        if lp
            if a_lt_b
                w.val  = log1p(-a / (a + b))
                w1.val = log(a / (a + b))
            else
                w.val = log(b / (a + b))
                w1.val = log1p(-b / (a + b))
            end
        else
            w.val = b / (a + b)
            w1.val = a / (a + b)
        end
        return
    end

    if min(a, b) <= 1.0
        do_swap = x > 0.5
        if do_swap
            a0 = b;  x0 = y
            b0 = a;  y0 = x
        else
            a0 = a;  x0 = x
            b0 = b;  y0 = y
        end

        if b0 < min(eps, eps * a0)
            w.val = fpser(a0, b0, x0, eps, lp)
            w1.val = lp ? Log1Exp(w.val) : 1.0 - w.val
            @goto L_end
        end

        if a0 < min(eps, eps * b0) && b0 * x0 <= 1.0
            w1.val = apser(a0, b0, x0, eps)
            @goto L_end_from_w1
        end

        did_bup = false
        if max(a0, b0) > 1.0
            b0 <= 1.0 && @goto L_w_bpser
            x0 >= 0.29 && @goto L_w1_bpser
            (x0 < 0.1 && (x0 * b0)^a0 <= 0.7) && @goto L_w_bpser

            if b0 > 15.0
                w1.val = 0.0
                @goto L131
            end
        else
            a0 >= min(0.2, b0) && @goto L_w_bpser
            x0 ^ a0 <= 0.9 && @goto L_w_bpser
            x0 >= 0.3 && @goto L_w1_bpser
        end

        n = 20
        w1.val = bup(b0, a0, y0, x0, n, eps, false)
        did_bup = true
        b0 += n

        @label L131
        bgrat(b0, a0, y0, x0, w1, 15.0 * eps, ierr1, false)

        if w1.val == 0 || (0 < w1.val && w1.val < typemin(Float64))
            if did_bup
                w1.val = bup(b0-n, a0, y0, x0, n, eps, true)
            else
                w1.val = -Inf
            end

            bgrat(b0, a0, y0, x0, w1, 15.0 * eps, ierr1, true)
            if ierr1.val != 0
                ierr.val = 10 + ierr1.val
            end
            @goto L_end_from_w1_log
        end

        if ierr1.val != 0
            ierr.val = 10 + ierr1.val
        end

        if w1 < 0 
            println("bratio(a = $(a), b = $(b), x = $(x)): bgrat() -> w1 = $(w1)")
        end
        @goto L_end_from_w1
    else
        λ = isfinite(a + b) ? ((a > b) ? (a + b) * y - b : a - (a + b) * x) : a * y - b * x
        do_swap = λ < 0.0
        if do_swap
            λ = -λ
            a0 = b;  x0 = y
            b0 = a;  y0 = x
        else
            a0 = a;  x0 = x
            b0 = b;  y0 = y
        end

        if b0 < 40.0
            if b0 * x0 <= 0.7 || (lp && λ > 650.0)
                @goto L_w_bpser
            else
                @goto L140
            end
        elseif a0 > b0
            if b0 <= 100.0 || λ > b0 * 0.03
                @goto L_bfrac
            end
        elseif a0 <= 100.0
            @goto L_bfrac
        elseif λ > a0 * 0.03
            @goto L_bfrac
        end

        w.val = basym(a0, b0, λ, eps * 100.0, lp)
        w1.val = lp ? Log1Exp(w.val) : 1.0 - w.val
        @goto L_end
    end

    @label L_w_bpser
    w.val = bpser(a0, b0, x0, eps, lp)
    w1.val = lp ? Log1Exp(w.val) : 1.0 - w.val
    @goto L_end

    @label L_w1_bpser
    w1.val = bpser(b0, a0, y0, eps, lp)
    w.val  = lp ? Log1Exp(w1.val) : 1.0 - w1.val
    @goto L_end

    @label L_bfrac
    w.val = bfrac(a0, b0, x0, y0, λ, eps * 15.0, lp)
    w1.val = lp ? Log1Exp(w.val) : 1.0 - w.val
    @goto L_end

    @label L140
    n = Int64(trunc(b0))
    b0 -= n
    if b0 == 0.0
        n -= 1
        b0 = 1.0
    end

    w.val = bup(b0, a0, y0, x0, n, eps, false)
    if w.val < typemin(Float64) && lp
        b0 += n
        @goto L_w_bpser
    end

    if x0 <= 0.7
        w.val += bpser(a0, b0, x0, eps, false)
        @goto L_end_from_w
    end

    if a0 <= 15.0
        n = 20
        w.val += bup(a0, b0, x0, y0, n, eps, false)
        a0 += n
    end

    bgrat(a0, b0, x0, y0, w, 15.0 * eps, ierr1, false)
    if ierr1.val != 0
        ierr.val = 10 + ierr1.val
    end
    @goto L_end_from_w

    # TERMINATION OF THE PROCEDURE

    @label L200
    if a == 0.0
        ierr.val = 6
        return
    end

    @label L201
    w.val = D0(lp); w1.val = D1(lp)
    return

    @label L210
    if b == 0.0
        ierr.val = 7
        return
    end

    @label L211
    w.val = D1(lp); w1.val = D0(lp)
    return
    
    @label L_end_from_w
    if lp
        w1.val = log1p(-w.val)
        w.val = log(w.val)
    else
        w1.val = 1.0 - w.val
    end
    @goto L_end

    @label L_end_from_w1
    if lp
        w.val = log1p(-w1.val)
        w1.val = log(w1.val)
    else
        w.val = 1.0 - w1.val
    end
    @goto L_end

    @label L_end_from_w1_log
    if lp
        w.val = Log1Exp(w1.val)
    else
        w.val = -expm1(w1.val)
        w1.val = exp(w1.val)
    end
    @goto L_end
    
    @label L_end
    if do_swap
        w.val, w1.val = w1.val, w.val
    end
    return
end # bratio

function fpser(a::Real, b::Real, x::Real, eps::Real, lp::Bool)
    ans, c, t = 0.0, 0.0, 0.0

    if lp
        ans = a * log(x)
    elseif a > eps * 0.001
        t = a * log(x)
        if t < exparg(1)
            return 0.0
        end
        ans = exp(t)
    else
        ans = 1.0
    end

    if lp
        ans += log(b) - log(a)
    else
        ans *= b / a
    end

    tol = eps / a
    an = a + 1.0
    t = x
    s = t / an

    while true
        an += 1.0
        t = x * t
        c = t / an
        s += c
        if abs(c) <= tol
            break
        end
    end

    if lp
        ans += log1p(a * s)
    else
        ans *= a * s + 1.0
    end
    return ans
end # fpser

function apser(a::Real, b::Real, x::Real, eps::Real)
    g = 0.577215664901533
    c, aj = 0.0, 0.0
    bx = b * x
    t = x - bx

    if b * eps <= 0.02
        c = log(x) + psi(b) + g + t
    else
        c = log(bx) + g + t
    end

    tol = eps * 5.0 * abs(c)
    j, s = 1.0, 0.0
    while true
        j += 1.0
        t *= x - bx / j
        aj = t / j
        s += aj
        if abs(aj) <= tol
            break
        end
    end
    return -a * (c + s)
end # apser

function bpser(a::Real, b::Real, x::Real, eps::Real, lp::Bool)
    i, m = 0, 0
    ans, c, t, u, z, b0, apb = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0

    if x == 0.0
        return D0(lp)
    end

    a0 = min(a, b)
    if a0 >= 1.0
        z = a * log(x) - betaln(a, b)
        ans = lp ? z - log(a) : exp(z) / a
    else
        b0 = max(a, b)

        if b0 < 8.0
            if b0 <= 1.0
                if lp
                    ans = a * log(x)
                else
                    ans = x ^ a
                    if ans == 0.0
                        return ans
                    end
                end

                apb = a + b
                if apb > 1.0
                    u = a + b - 1.0
                    z = (gam1(u) + 1.0) / apb
                else
                    z = gam1(apb) + 1.0
                end

                c = (gam1(a) + 1.0) * (gam1(b) + 1.0) / z
                if lp
                    ans += log(c * (b / apb))
                else
                    ans *=  c * (b / apb)
                end
            else
                u = gamln1(a0)
                m = Int64(trunc(b0 - 1.0))
                if m >= 1
                    c = 1.0
                    for i = 1:m
                        b0 -= 1.0
                        c *= b0 / (a0 + b0)
                    end
                    u += log(c)
                end

                z = a * log(x) - u
		        b0 -= 1.0
                apb = a0 + b0
                if apb > 1.0
                    u = a0 + b0 - 1.0
                    t = (gam1(u) + 1.0) / apb
                else
                    t = gam1(apb) + 1.0
                end

                if lp
                    ans = z + log(a0 / a) + log1p(gam1(b0)) - log(t)
                else
                    ans = exp(z) * (a0 / a) * (gam1(b0) + 1.0) / t
                end
            end
        else
            u = gamln1(a0) + algdiv(a0, b0)
            z = a * log(x) - u
            
            if lp
                ans = z + log(a0 / a)
            else
                ans = a0 / a * exp(z)
            end
        end
    end

    if ans == D0(lp) || (!lp && a <= eps * 0.1)
        return ans
    end

    tol = eps / a
    n, sum, w, c = 0.0, 0.0, 0.0, 1.0
    while true
        n += 1.0
        c *= (1.0 - b / n) * x
        w = c / (a + n)
        sum += w
        if n >= 1e7 || abs(w) <= tol
            break
        end
    end

    if abs(w) > tol
        tmp1 = (lp && !(a*sum > -1.0 && abs(log1p(a * sum)) < eps * abs(ans)))
        tmp2 = (!lp && abs(a*sum + 1.0) != 1.0)
        if tmp1 || tmp2
            print("bpser(a = $(a), b = $(b), x = $(x),...) did not converge ")
            println("(n = 1e7, |w|/tol = $(abs(w)/tol) > 1; A = $(ans)")
        end
    end
        
    if lp
        if a * sum > -1.0
            ans += log1p(a * sum)
        else
            if ans > -Inf
                print("pbeta(*, log_p = true) -> bpser(a = $(a), b = $(b), x = $(x),...)")
                println(" underflow to -Inf")
            end
            ans = -Inf
        end
    elseif a * sum > -1.0
        ans *= (a * sum + 1.0)
    else
        ans = 0.0
    end
    return ans
end # bpser

function bup(a::Real, b::Real, x::Real, y::Real, n::Int64, eps::Real, lp::Bool)
    i, k, mu = 0, 0, 0
    ret_val, d, l = 0.0, 0.0, 0.0

    apb = a + b
    ap1 = a + 1.0
    if n > 1 && a >= 1. && apb >= ap1 * 1.1
        mu = Int64(trunc(abs(exparg(1))))
        k = Int64(trunc(exparg(0)))
        if mu > k
            mu = k
        end
        d = exp(-mu)
    else
        mu = 0
        d = 1.0
    end

    ret_val = lp ? brcmp1(mu, a, b, x, y, true) - log(a) : brcmp1(mu, a, b, x, y, false)  / a
    if n == 1 || (lp && ret_val == -Inf) || (!lp && ret_val == 0.0)
        return ret_val
    end

    nm1 = n - 1
    w = d
    k = 0
    if b > 1.0
        if y > 1e-4
            r = (b - 1.0) * x / y - a
            if r >= 1.0
                k = r < nm1 ? Int64(trunc(r)) : nm1
            end
        else
            k = nm1
        end

        for i = 0:(k-1)
            l = Float64(i)
            d *= (apb + l) / (ap1 + l) * x
            w += d
        end
    end

    for i = k:(nm1-1)
        l = Float64(i)
        d *= (apb + l) / (ap1 + l) * x
        w += d
        if d <= eps * w
            break
        end
    end
        
    if lp
        ret_val += log(w)
    else
        ret_val *= w
    end
    return ret_val
end # bup

function bfrac(a::Real, b::Real, x::Real, y::Real, lambda::Real, eps::Real, lp::Bool)
    λ = lambda
    e, t, w, r0, beta, alpha = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    
    if isinf(λ)
        return NaN
    end

    brc = brcomp(a, b, x, y, lp)
    if isnan(brc)
        return NaN
    end
    if !lp && brc == 0.0
        return 0.0
    end

    c = λ + 1.0
    c0 = b / a
    c1 = 1.0 / a + 1.0
    yp1 = y + 1.0

    n = 0.0
    p = 1.0
    s = a + 1.0
    an = 0.0
    bn = 1.0
    anp1 = 1.0
    bnp1 = c / c1
    r = c1 / c

    while true
        n += 1.0
        t = n / a
        w = n * (b - n) * x
        e = a / s
        alpha = p * (p + c0) * e * e * (w * x)
        e = (t + 1.0) / (c1 + t + t)
        beta = n + w / s + e * (c + n * yp1)
        p = t + 1.0
        s += 2.0

        t = alpha * an + beta * anp1;	an = anp1;	anp1 = t
        t = alpha * bn + beta * bnp1;	bn = bnp1;	bnp1 = t

        r0 = r
        r = anp1 / bnp1

        if abs(r - r0) <= eps * r
            break
        end

        an /= bnp1
        bn /= bnp1
        anp1 = r
        bnp1 = 1.0
        if n >= 10000
            break
        end
    end

    if n >= 10000 && abs(r - r0) > eps * r
        print("bfrac(a = $(a), b = $(b), x = $(x), y = $(y), λ = $(λ))")
        println(" did *not* converge (in 10000 steps)")
    end
    return lp ? brc + log(r) : brc * r
end # bfrac

function brcomp(a::Real, b::Real, x::Real, y::Real, lp::Bool)
    const__ = 0.398942280401433
    i, n = 0, 0
    c, e, u, v, z, b0, apb = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0

    if x == 0.0 || y == 0.0
        return D0(lp)
    end

    a0 = min(a, b)
    if a0 < 8.0
        lnx, lny = 0.0
        if x <= 0.375
            lnx = log(x)
            lny = alnrel(-x)
        else
            if y > 0.375
                lnx = log(x)
                lny = log(y)
            else
                lnx = alnrel(-y)
                lny = log(y)
            end
        end

        z = a * lnx + b * lny
        if a0 >= 1.0
            z -= betaln(a, b)
            return Dexp(z, lp)
        end

        b0 = max(a, b)
        if b0 >= 8.0
            u = gamln1(a0) + algdiv(a0, b0)
            return lp ? log(a0) + (z - u)  : a0 * exp(z - u)
        end

        if b0 <= 1.0
            e_z = Dexp(z, lp)
            if !lp && e_z == 0.0
                return 0.0
            end

            apb = a + b
            if apb > 1.0
                u = a + b - 1.0
                z = (gam1(u) + 1.0) / apb
            else
                z = gam1(apb) + 1.0
            end

            c = (gam1(a) + 1.) * (gam1(b) + 1.) / z
            return lp ? e_z + log(a0 * c) - log1p(a0 / b0) : e_z * (a0 * c) / (a0 / b0 + 1.0)
        end

        u = gamln1(a0)
        n = Int64(trunc((b0 - 1.0)))
        if n >= 1
            c = 1.0
            for i = 1:n
                b0 -= 1.0
                c *= b0 / (a0 + b0)
            end
            u = log(c) + u
        end

        z -= u
        b0 -= 1.0
        apb = a0 + b0
        t = 0.0
        if apb > 1.0
            u = a0 + b0 - 1.0
            t = (gam1(u) + 1.0) / apb
        else
            t = gam1(apb) + 1.0
        end

        return lp ? log(a0) + z + log1p(gam1(b0))  - log(t) : a0 * exp(z) * (gam1(b0) + 1.0) / t
    else
        h, x0, y0, λ = 0.0, 0.0, 0.0, 0.0
        if a <= b
            h = a / b;
            x0 = h / (h + 1.0)
            y0 = 1.0 / (h + 1.0)
            λ = a - (a + b) * x
        else
            h = b / a;
            x0 = 1.0 / (h + 1.0)
            y0 = h / (h + 1.0)
            λ = (a + b) * y - b
        end

        e = -λ / a
        if abs(e) > 0.6
            u = e - log(x / x0)
        else
            u = rlog1(e)
        end

        e = λ / b
        if abs(e) <= 0.6
            v = rlog1(e)
        else
            v = e - log(y / y0)
        end

        z = lp ? -(a * u + b * v) : exp(-(a * u + b * v))
        return lp ? -log(√(2.0 * π)) + 0.5 * log(b * x0) + z - bcorr(a, b) : const__ * √(b * x0) * z * exp(-bcorr(a, b))
    end
end # brcomp

function brcmp1(mu::Int64, a::Real, b::Real, x::Real, y::Real, lp::Bool)
    const__ = 0.398942280401433
    c, t, u, v, z, b0, apb = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0

    a0 = min(a, b)
    if a0 < 8.0
        lnx, lny = 0.0, 0.0
        if x <= 0.375
            lnx = log(x)
            lny = alnrel(-x)
        elseif y > 0.375
            lnx = log(x)
            lny = log(y)
        else
            lnx = alnrel(-y)
            lny = log(y)
        end

        z = a * lnx + b * lny
        if a0 >= 1.0
            z -= betaln(a, b)
            return esum(mu, z, lp)
        end

        b0 = max(a, b)
        if b0 >= 8.0
            u = gamln1(a0) + algdiv(a0, b0)
            return lp ? log(a0) + esum(mu, z - u, true) : a0  * esum(mu, z - u, false)
        elseif b0 <= 1.0
            ans = esum(mu, z, lp)
            if ans == lp ? -Inf : 0.0
                return ans
            end

            apb = a + b
            if apb > 1.0
                u = a + b - 1.0
                z = (gam1(u) + 1.0) / apb
            else
                z = gam1(apb) + 1.0
            end

            c = lp ? log1p(gam1(a)) + log1p(gam1(b)) - log(z) : (gam1(a) + 1.0) * (gam1(b) + 1.0) / z
            return lp ? ans + log(a0) + c - log1p(a0 / b0) : ans * (a0 * c) / (a0 / b0 + 1.0)
        end

        u = gamln1(a0)
        n = Int64(trunc(b0 - 1.0))
        if n >= 1.0
            c = 1.0
            for i = 1:n
                b0 -= 1.0
                c *= b0 / (a0 + b0)
            end
            u += log(c)
        end

        z -= u
	    b0 -= 1.0
        apb = a0 + b0
        if apb > 1.0
            t = (gam1(apb - 1.0) + 1.0) / apb
        else
            t = gam1(apb) + 1.0
        end

        return lp ? log(a0)+ esum(mu, z, true) + log1p(gam1(b0)) - log(t) : a0 * esum(mu, z, false) * (gam1(b0) + 1.0) / t
    else
        h, x0, y0, λ = 0.0, 0.0, 0.0, 0.0
        if a > b
            h = b / a
            x0 = 1. / (h + 1.0)
            y0 = h / (h + 1.0)
            λ = (a + b) * y - b
        else
            h = a / b
            x0 = h / (h + 1.0)
            y0 = 1. / (h + 1.0)
            λ = a - (a + b) * x
        end

        lx0 = -log1p(b / a)
        e = -λ / a
        if abs(e) > 0.6
            u = e - log(x / x0)
        else
            u = rlog1(e)
        end

        e = λ / b
        if abs(e) > 0.6
            v = e - log(y / y0)
        else
            v = rlog1(e)
        end

        z = esum(mu, -(a * u + b * v), lp)
        return lp ? log(const__) + (log(b) + lx0) / 2.0 + z - bcorr(a, b) : const__ * √(b * x0) * z * exp(-bcorr(a, b))
    end
end # brcmp1

function bgrat(a::Real, b::Real, x::Real, y::Real, w::RefFloat, eps::Real, ierr::RefInt, log_w::Bool)
    n_terms_bgrat = 30
    bm1 = b - 1.0
    nu = a + bm1 * 0.5
    lnx = (y > 0.375) ? log(x) : alnrel(-y)
    z = -nu * lnx
    c = [0.0 for _ = 1:n_terms_bgrat]
    d = [0.0 for _ = 1:n_terms_bgrat]

    if b * z == 0.0
        print("bgrat(a = $(a), b = $(b), x = $(x), y = $(y)): z = $(z) ")
        println("b * z == 0 underflow, hence inacurrate pbeta()")
        ierr.val = 1
        return
    end

    log_r = log(b) + log1p(gam1(b)) + b * log(z) + nu * lnx
    log_u = log_r - (algdiv(b, a) + b * log(nu))
    u = exp(log_u)

    if log_u == -Inf
        ierr.val = 2
        return
    end

    u_0 = u == 0.0
    l = log_w ? ((w == -Inf) ? 0.0 : exp(w - log_u)) : ((w == 0.0) ? 0.0 : exp(log(w) - log_u))
    q_r = grat_r(b, z, log_r, eps)
    v = 0.25 / (nu * nu)
    t2 = lnx * 0.25 * lnx
    j = q_r
    sum = j
    t, cn, n2 = 1.0, 1.0, 0.0
    for n = 1:n_terms_bgrat
        bp2n = b + n2
        j = (bp2n * (bp2n + 1.0) * j + (z + bp2n + 1.0) * t) * v
		n2 += 2.0
		t *= t2
        cn /= n2 * (n2 + 1.0)
        nm1 = n - 1
        c[n] = cn
        s = 0.0
        if n > 1
            coef = b - n
            for i = 1:nm1
                s += coef * c[i] * d[n - i];
                coef += b;
            end
        end

        d[n] = bm1 * cn + s / n
        dj = d[n] * j
        sum += dj
        if sum <= 0.0
            ierr.val = 3
            return
        end

        if abs(dj) <= eps * (sum + l)
            ierr.val = 0
            break
        elseif n == n_terms_bgrat
            ierr.val = 4
            print("bgrat(a = $(a), b = $(b), x = $(x)) no convergence: ")
            println("NOTIFY! dj = $(dj), rel_err = $(abs(dj)/(sum + l))")
        end
    end

    if log_w
        w.val = logspace_add(w.val, log_u + log(sum))
    else
        w.val += (u_0 ? exp(log_u + log(sum)) : u * sum)
    end
    return
end # bgrat

function grat_r(a::Real, x::Real, log_r::Real, eps::Real)
    if  a * x == 0.0
        if x <= a
            return exp(-log_r)
        else
            return 0.0
        end
    elseif a == 0.5
        if x < 0.25
            p = erf__(√x)
            return (1.0 - p) * exp(-log_r)
        else
            sx = √x
            return erfc1(1, sx) / sx * √π
        end
    elseif x < 1.1
        an = 3.0
        c = x
        sum = x / (a + 3.0)
        tol = eps * 0.1 / (a + 1.0)
        t = 0

        while true
            an += 1.0
			c *= -(x / an)
			t = c / (a + an)
            sum += t
            if abs(t) <= tol
                break
            end
        end

        j = a * x * ((sum / 6.0 - 0.5 / (a + 2.0)) * x + 1.0 / (a + 1.0))
        z = a * log(x)
        h = gam1(a)
        g = h + 1.0

        if (x >= 0.25 && (a < x / 2.59)) || (z > -0.13394)
            l = rexpm1(z)
            q = ((l + 1.0) * j - l) * g - h
            if q <= 0.0
                return 0.0
            else
                return q * exp(-log_r)
            end
        else
            p = exp(z) * g * (1.0 - j)
            return (1.0 - p) * exp(-log_r)
        end
    else
        a2n_1, a2n, b2n_1, b2n = 1.0, 1.0, x, x + (1.0 - a)
        c, am0, an0 = 1.0, 0.0, 0.0

        while true
            a2n_1 = x * a2n + c * a2n_1
			b2n_1 = x * b2n + c * b2n_1
			am0 = a2n_1 / b2n_1
			c += 1.0
			c_a = c - a
			a2n = a2n_1 + c_a * a2n
			b2n = b2n_1 + c_a * b2n
            an0 = a2n / b2n
            
            if abs(an0 - am0) < eps * an0
                break
            end
        end
        return an0
    end
end # grat_r

function basym(a::Real, b::Real, lambda::Real, eps::Real, lp::Bool)
    λ = lambda
    num_IT = 20
    e0 = 1.12837916709551
    e1 = 0.353553390593274
    ln_e0 = 0.120782237635245

    a0 = [0.0 for _ = 1:(num_IT+1)]
    b0 = [0.0 for _ = 1:(num_IT+1)]
    c = [0.0 for _ = 1:(num_IT+1)]
    d = [0.0 for _ = 1:(num_IT+1)]

    f = a * rlog1(-λ/a) + b * rlog1(λ/b)
    t = 0.0

    if lp
        t = -f
    else
        t = exp(-f)
        if t == 0.0
            return 0.0
        end
    end

    z0 = √f
	z = z0 / e1 * 0.5
	z2 = f + f
    h, r0, r1, w0 = 0.0, 0.0, 0.0, 0.0
    
    if a < b
        h = a / b
        r0 = 1.0 / (h + 1.0)
        r1 = (b - a) / b
        w0 = 1.0 / √(a * (h + 1.0))
    else
        h = b / a
        r0 = 1.0 / (h + 1.0)
        r1 = (b - a) / a
        w0 = 1.0 / √(b * (h + 1.0))
    end

    a0[1] = r1 * 0.66666666666666663
    c[1] = a0[1] * -0.5
    d[1] = -c[1]
    j0 = 0.5 / e0 * erfc1(1, z0)
    j1 = e1
    sum = j0 + d[1] * w0 * j1

    s, h2, hn, w, znm1, zn = 1.0, h * h, 1.0, w0, z, z2
    for n = 2:num_IT
        hn *= h2
        a0[n] = r0 * 2.0 * (h * hn + 1.0) / (n + 2.0)
        np1 = n + 1
        s += hn
        a0[np1] = r1 * 2.0 * s / (n + 3.0)

        for i = n:np1
            r = (i + 1.0) * -0.5
            b0[1] = r * a0[1]
            for m = 2:i
                bsum = 0.0
                for j = 1:(m-1)
                    mmj = m - j
                    bsum += (j * r - mmj) * a0[j] * b0[mmj]
                end
                b0[m] = r * a0[m] + bsum / m
            end
            
            c[i] = b0[i] / (i + 1.0)
            dsum = 0.0
            for j = 1:(i-1)
                dsum += d[i - j] * c[j]
            end

            d[i] = -(dsum + c[i])
        end

        j0 = e1 * znm1 + (n - 1.0) * j0
		j1 = e1 * zn + n * j1
		znm1 = z2 * znm1
		zn = z2 * zn
		w *= w0
		t0 = d[n] * w * j0
		w *= w0
		t1 = d[np1] * w * j1
        sum += t0 + t1
        if abs(t0) + abs(t1) <= eps * sum
            break
        end
    end

    if lp
        return ln_e0 + t - bcorr(a, b) + log(sum)
    else
        u = exp(-bcorr(a, b))
        return e0 * t * u * sum
    end
end # basym

function exparg(l::Int64)
    lnb = 0.69314718055995
    m = l == 0 ? 1024 : -1022
    return m * lnb * 0.99999
end

function esum(mu::Int64, x::Real, lp::Bool)
    if lp
        return x + mu
    end

    w = 0.0
    ret = exp(mu) * exp(x)
    if x > 0.0
        if mu > 0
            return ret
        end
        w = mu + x
        if w < 0.0
            return ret
        end
    else
        if mu < 0
            return ret
        end
        w = mu + x
        if w > 0.0
            return ret
        end
    end
    return exp(w)
end # esum

function rexpm1(x::Real)
    p1 = 9.14041914819518e-10
    p2 = 0.0238082361044469
    q1 = -0.499999999085958
    q2 = 0.107141568980644
    q3 = -0.0119041179760821
    q4 = 5.95130811860248e-4

    if abs(x) <= 0.15
        return x * (((p2 * x + p1) * x + 1.0) / ((((q4 * x + q3) * x + q2) * x + q1) * x + 1.0))
    else
        w = exp(x)
        if x > 0.0
            return w * (1.0 - 1.0 / w)
        else
            return return w - 1.0
        end
    end
end # rexpm1

function alnrel(a::Real)
    if abs(a) > 0.375
        return log(1.0 + a)
    end

    p1 = -1.29418923021993
	p2 = 0.405303492862024
	p3 = -0.0178874546012214
	q1 = -1.62752256355323
	q2 = 0.747811014037616
    q3 = -0.0845104217945565
    
    t = a / (a + 2.)
	t2 = t * t
    w = (((p3 * t2 + p2) * t2 + p1) * t2 + 1.0) / (((q3 * t2 + q2) * t2 + q1) * t2 + 1.0)
    return t * 2.0 * w
end # alnrel

function rlog1(x::Real)
    a = 0.0566749439387324
    b = 0.0456512608815524
    p0 = 0.333333333333333
    p1 = -0.224696413112536
    p2 = 0.00620886815375787
    q1 = -1.27408923933623
    q2 = 0.354508718369557

    h, r, t, w, w1 = 0.0, 0.0, 0.0, 0.0, 0.0
    if x < -0.39 || x > 0.57
        w = x + 1.0
        return x - log(w)
    end

    if x < -0.18
        h = x + 0.3
        h /= 0.7
        w1 = a - h * 0.3
    elseif x > 0.18
        h = x * 0.75 - 0.25;
        w1 = b + h / 3.0
    else
        h = x
        w1 = 0.0
    end

    r = h / (h + 2.0)
    t = r * r
    w = ((p2 * t + p1) * t + p0) / ((q2 * t + q1) * t + 1.0)
    return t * 2.0 * (1.0 / (1.0 - r) - r * w) + w1
end

function erf__(x::Real)
    c = 0.564189583547756
    a = [7.7105849500132e-5, -0.00133733772997339,
        0.0323076579225834, 0.0479137145607681,
        0.128379167095513]
    b = [0.00301048631703895, 0.0538971687740286,
        0.375795757275549 ]
    p = [-1.36864857382717e-7, 0.564195517478974,
        7.21175825088309, 43.1622272220567, 152.98928504694,
        339.320816734344, 451.918953711873, 300.459261020162]
    q = [1.0, 12.7827273196294, 77.0001529352295,
        277.585444743988, 638.980264465631, 931.35409485061,
        790.950925327898, 300.459260956983]
    r = [2.10144126479064, 26.2370141675169,
        21.3688200555087, 4.6580782871847, 0.282094791773523]
    s = [94.153775055546, 187.11481179959,
        99.0191814623914, 18.0124575948747]
    
    t, x2, bot, top = 0.0, 0.0, 0.0, 0.0
    ax = abs(x)
    if ax <= 0.5
        t = x * x
        top = (((a[1] * t + a[2]) * t + a[3]) * t + a[4]) * t + a[5] + 1.0
        bot = ((b[1] * t + b[2]) * t + b[3]) * t + 1.0
        return x * (top / bot)
    end

    if ax <= 4.0
        top = ((((((p[1] * ax + p[2]) * ax + p[3]) * ax + p[4]) * ax + p[5]) * ax
		    + p[6]) * ax + p[7]) * ax + p[8]
        bot = ((((((q[1] * ax + q[2]) * ax + q[3]) * ax + q[4]) * ax + q[5]) * ax
            + q[6]) * ax + q[7]) * ax + q[8]
        ret = 0.5 - exp(-x * x) * top / bot + 0.5
        return (x < 0) ? -ret : ret
    end

    if ax >= 5.8
        return x > 0 ? 1 : -1
    end

    x2 = x * x
    t = 1.0 / x2
    top = (((r[1] * t + r[2]) * t + r[3]) * t + r[4]) * t + r[5]
    bot = (((s[1] * t + s[2]) * t + s[3]) * t + s[4]) * t + 1.0
    t = (c - top / (x2 * bot)) / ax
    ret = 1.0 - exp(-x2) * t
    return (x < 0) ? -ret : ret
end # erf__

function erfc1(ind::Int64, x::Real)
    c = 0.564189583547756
    a = [7.7105849500132e-5, -0.00133733772997339,
        0.0323076579225834, 0.0479137145607681,
        0.128379167095513]
    b = [0.00301048631703895, 0.0538971687740286,
        0.375795757275549 ]
    p = [-1.36864857382717e-7, 0.564195517478974,
        7.21175825088309, 43.1622272220567, 152.98928504694,
        339.320816734344, 451.918953711873, 300.459261020162]
    q = [1.0, 12.7827273196294, 77.0001529352295,
        277.585444743988, 638.980264465631, 931.35409485061,
        790.950925327898, 300.459260956983]
    r = [2.10144126479064, 26.2370141675169,
        21.3688200555087, 4.6580782871847, 0.282094791773523]
    s = [94.153775055546, 187.11481179959,
        99.0191814623914, 18.0124575948747]

    ret_val, e, t, w, bot, top = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    ax = abs(x)

    if ax <= 0.5
        t = x * x
		top = (((a[1] * t + a[2]) * t + a[3]) * t + a[4]) * t + a[5] + 1.0
		bot = ((b[1] * t + b[2]) * t + b[3]) * t + 1.0
        ret_val = 0.5 - x * (top / bot) + 0.5
        if ind != 0
            ret_val = exp(t) * ret_val
        end
        return ret_val
    end

    if ax <= 4.0
        top = ((((((p[1] * ax + p[2]) * ax + p[3]) * ax + p[4]) * ax + p[5]) * ax
			+ p[6]) * ax + p[7]) * ax + p[8]
		bot = ((((((q[1] * ax + q[2]) * ax + q[3]) * ax + q[4]) * ax + q[5]) * ax
			+ q[6]) * ax + q[7]) * ax + q[8]
        ret_val = top / bot
    else
        if x <= -5.6
            ret_val = 2.0
            if ind != 0
                ret_val = exp(x * x) * 2.0
            end
            return ret_val
        end

        if ind == 0 && (x > 100.0 || x * x > -exparg(1))
            return 0.0
        end

        t = 1.0 / (x * x)
		top = (((r[1] * t + r[2]) * t + r[3]) * t + r[4]) * t + r[5]
		bot = (((s[1] * t + s[2]) * t + s[3]) * t + s[4]) * t + 1.0
        ret_val = (c - t * top / bot) / ax
    end

    if ind != 0
        if x < 0
            ret_val = exp(x * x) * 2.0 - ret_val
        end
    else
        w = x * x
		t = w
		e = w - t
        ret_val = (1.0 - e) * exp(-t) * ret_val
        if x < 0.0
            ret_val = 2.0 - ret_val
        end
    end
    return ret_val
end # erfc1

function gam1(a::Real)
    w, bot, top = 0.0, 0.0, 0.0
    t = a
    d = a - 0.5

    if d > 0.0
        t = d - 0.5
    end

    if t < 0.0
        r = [-0.422784335098468, -0.771330383816272,
            -0.244757765222226, 0.118378989872749, 9.30357293360349e-4,
            -0.0118290993445146, 0.00223047661158249, 2.66505979058923e-4,
            -1.32674909766242e-4]
        s1 = 0.273076135303957
        s2 = 0.0559398236957378
        top = (((((((r[9] * t + r[8]) * t + r[7]) * t + r[6]) * t + r[5]) * t + r[4]) * t + r[3]) * t + r[2]) * t + r[1]
        bot = (s2 * t + s1) * t + 1.0
        w = top / bot
        if d > 0.0
            return t * w / a
        else
            return a * (w + 1.0)
        end
    elseif t == 0
        return 0.0
    else
        p = [0.577215664901533, -0.409078193005776,
            -0.230975380857675, 0.0597275330452234, 0.0076696818164949,
            -0.00514889771323592, 5.89597428611429e-4]
        q = [1.0, 0.427569613095214, 0.158451672430138,
            0.0261132021441447, 0.00423244297896961]

        top = (((((p[7] * t + p[6]) * t + p[5]) * t + p[4]) * t + p[3]) * t + p[2]) * t + p[1]
        bot = (((q[5] * t + q[4]) * t + q[3]) * t + q[2]) * t + 1.0
        w = top / bot
        
        if d > 0.0
            return t / a * (w - 1.0)
        else
            return a * w
        end
    end
end # gam1

function gamln1(a::Real)
    w = 0.0
    if a < 0.6
        p0 = 0.577215664901533
        p1 = 0.844203922187225
        p2 = -0.168860593646662
        p3 = -0.780427615533591
        p4 = -0.402055799310489
        p5 = -0.0673562214325671
        p6 = -0.00271935708322958
        q1 = 2.88743195473681
        q2 = 3.12755088914843
        q3 = 1.56875193295039
        q4 = 0.361951990101499
        q5 = 0.0325038868253937
        q6 = 6.67465618796164e-4
        w = ((((((p6 * a + p5)* a + p4)* a + p3)* a + p2)* a + p1)* a + p0) /
            ((((((q6 * a + q5)* a + q4)* a + q3)* a + q2)* a + q1)* a + 1.0)
        return -a * w
    else
        r0 = 0.422784335098467
        r1 = 0.848044614534529
        r2 = 0.565221050691933
        r3 = 0.156513060486551
        r4 = 0.017050248402265
        r5 = 4.97958207639485e-4
        s1 = 1.24313399877507
        s2 = 0.548042109832463
        s3 = 0.10155218743983
        s4 = 0.00713309612391
        s5 = 1.16165475989616e-4
        x = a - 0.5 - 0.5
        w = (((((r5 * x + r4) * x + r3) * x + r2) * x + r1) * x + r0) /
            (((((s5 * x + s4) * x + s3) * x + s2) * x + s1) * x + 1.0)
        return x * w
    end
end # gamln1

function psi(x::Real)
    piov4 = 0.785398163397448
    dx0 = 1.461632144968362341262659542325721325
    p1 = [0.0089538502298197, 4.77762828042627,
        142.441585084029, 1186.45200713425, 3633.51846806499,
        4138.10161269013, 1305.60269827897]
    q1 = [44.8452573429826, 520.752771467162,
        2210.0079924783, 3641.27349079381, 1908.310765963,
        6.91091682714533e-6]
    p2 = [-2.12940445131011, -7.01677227766759,
        -4.48616543918019, -0.648157123766197]
    q2 = [32.2703493791143, 89.2920700481861,
        54.6117738103215, 7.77788548522962]

    i, m, n, nq = 0, 0, 0, 0
    w, z = 0.0, 0.0, 0.0
    den, sgn, xmx0, upper = 0.0, 0.0, 0.0, 0.0

    xmax1 = Float64(typemax(Int32))
    d2 = 1.0 / eps()
    if xmax1 > d2
        xmax1 = d2
    end

    xsmall = 1e-9
    aug = 0.0
    if x < 0.5
        if abs(x) <= xsmall
            if x == 0.0
                return 0.0
            end
            aug = -1.0 / x
        else
            w = -x
            sgn = piov4
            if w <= 0.0
                w = -w
                sgn = -sgn
            end

            if w >= xmax1
                return 0.0
            end

            nq =  Int64(trunc(w))
            w -= Float64(nq)
            nq = Int64(trunc(w * 4.0))
            w = (w - Float64(nq) * 0.25) * 4.0
            n = nq / 2
            if n + n != nq
                w = 1.0 - w
            end

            z = piov4 * w
            m = n / 2
            if m + m != n
                sgn = -sgn
            end

            n = (nq + 1) / 2
            m = n / 2
            m += m
            if m == n
                if z == 0.0
                    return 0.0
                end
                aug = sgn * (cos(z) / sin(z) * 4.0)
            else
                aug = sgn * (sin(z) / cos(z) * 4.0)
            end
        end

        x = 1.0 - x
    end

    if x <= 3.0
        den = x
        upper = p1[1] * x
        for i = 1:5
            den = (den + q1[i]) * x
            upper = (upper + p1[i + 1]) * x
        end

        den = (upper + p1[7]) / (den + q1[6])
        xmx0 = x - dx0
        return den * xmx0 + aug
    end

    if x < xmax1
        w = 1.0 / (x * x)
		den = w
		upper = p2[1] * w

		for i = 1:3
			den = (den + q2[i]) * w
			upper = (upper + p2[i+1]) * w
        end
        aug = upper / (den + q2[4]) - 0.5 / x + aug
    end

    return aug + log(x)
end # psi

function betaln(a0::Real, b0::Real)
    e = 0.918938533204673
    a = min(a0, b0)
    b = max(a0, b0)
    n = 0

    if a < 8.0
        if a < 1.0
            if b < 8.0
                return gamln(a) + (gamln(b) - gamln(a+b))
            else
                return gamln(a) + algdiv(a, b)
            end
        end

        w = 0.0
        if a < 2.0
            if b <= 2.0
                return gamln(a) + gamln(b) - gsumln(a, b)
            end

            if b < 8.0
                w = 0.0
                @goto L40
            end
            return gamln(a) + algdiv(a, b)
        end

        if b <= 1e3
            n = Int64(trunc(a - 1.0))
            w = 1.0
            for i = 1:n
                a -= 1.0
                h = a / b
                w *= h / (h + 1.0)
            end
            w = log(w)

            if b >= 8.0
                return w + gamln(a) + algdiv(a, b)
            end

            @label L40
            n = Int64(trunc(b - 1.0))
            z = 1.0
            for i = 1:n
                b -= 1.0
                z *= b / (a + b)
            end
            return w + log(z) + (gamln(a) + (gamln(b) - gsumln(a, b)))
        else
            n = Int64(trunc(a - 1.0))
            w = 1.0
            for i = 1:n
                a -= 1.0
                w *= a / (a / b + 1.0)
            end
            return log(w) - n * log(b) + (gamln(a) + algdiv(a, b))
        end
    else
        w = bcorr(a, b)
	    h = a / b
	    u = -(a - 0.5) * log(h / (h + 1.0))
        v = b * alnrel(h)
        
        return log(b) * -0.5 + e + w - v - u
    end
end # betaln

function gsumln(a::Real, b::Real)
    x = a + b - 2.0
    if x <= 0.25
        return gamln1(x + 1.0)
    end

    if x <= 1.25
        return gamln1(x) + alnrel(x)
    end

    return gamln1(x - 1.0) + log(x * (x + 1.0))
end # gsumln

function bcorr(a0::Real, b0::Real)
    c0 = 0.0833333333333333
    c1 = -0.00277777777760991
    c2 = 7.9365066682539e-4
    c3 = -5.9520293135187e-4
    c4 = 8.37308034031215e-4
    c5 = -0.00165322962780713

    a = min(a0, b0)
    b = max(a0, b0)

    h = a / b
    c = h / (h + 1.0)
    x = 1.0 / (h + 1.0)
    x2 = x * x

    s3 = x + x2 + 1.0
    s5 = x + x2 * s3 + 1.0
    s7 = x + x2 * s5 + 1.0
    s9 = x + x2 * s7 + 1.0
    s11 = x + x2 * s9 + 1.0

    r1 = 1.0 / b
    t = r1 * r1
    w = ((((c5 * s11 * t + c4 * s9) * t + c3 * s7) * t + c2 * s5) * t + c1 * s3) * t + c0
    w *= c / b

    r1 = 1.0 / a
    t = r1 * r1
    ret_val = (((((c5 * t + c4) * t + c3) * t + c2) * t + c1) * t + c0) / a + w

    return ret_val
end # bcorr

function algdiv(a::Real, b::Real)
    c0 = 0.0833333333333333
    c1 = -0.00277777777760991
    c2 = 7.9365066682539e-4
    c3 = -5.9520293135187e-4
    c4 = 8.37308034031215e-4
    c5 = -0.00165322962780713

    h, c, x, d = 0.0, 0.0, 0.0, 0.0
    if a > b
        h = b / a
        c = 1.0 / (h + 1.0)
        x = h / (h + 1.0)
        d = a + (b - 0.5)
    else
        h = a / b
        c = h / (h + 1.0)
        x = 1.0 / (h + 1.0)
        d = b + (a - 0.5)
    end

    x2 = x * x
    s3 = x + x2 + 1.0
    s5 = x + x2 * s3 + 1.0
    s7 = x + x2 * s5 + 1.0
    s9 = x + x2 * s7 + 1.0
    s11 = x + x2 * s9 + 1.0

    t = 1.0 / (b * b)
    w = ((((c5 * s11 * t + c4 * s9) * t + c3 * s7) * t + c2 * s5) * t + c1 * s3) * t + c0
    w *= c / b

    u = d * alnrel(a / b)
    v = a * (log(b) - 1.0)

    return w - v - u
end # algdiv

function gamln(a::Real)
    d = 0.418938533204673
    c0 = 0.0833333333333333
    c1 = -0.00277777777760991
    c2 = 7.9365066682539e-4
    c3 = -5.9520293135187e-4
    c4 = 8.37308034031215e-4
    c5 = -0.00165322962780713

    if a <= 0.8
        return gamln1(a) - log(a)
    elseif a <= 2.25
        return gamln1(a - 1.0)
    elseif a < 10.0
        n = Int64(trunc(a - 1.25))
        t = a
        w = 1.0
        for i = 1:n
            t -= 1.0
            w *= t
        end
        return gamln1(t - 1.0) + log(w)
    else
        t = 1.0 / (a * a)
        w = (((((c5 * t + c4) * t + c3) * t + c2) * t + c1) * t + c0) / a
        return d + w + (a - 0.5) * (log(a) - 1.0)
    end
end # gamln
