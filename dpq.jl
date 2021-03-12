# lp: log_p, lt: lower_tail

function D0(lp)
    lp ? -Inf : 0.0
end

function D1(lp)
   lp ? 0.0 : 1.0 
end

function DT0(lt, lp)
    lt ? D0(lp) : D1(lp)
end

function DT1(lt, lp)
    lt ? D1(lp) : D0(lp)
end

function Dhalf(lp)
    lp ? -log(2.0) : 0.5
end

function DLval(p, lt)
    lt ? p : 1.0 - p;
end

function DCval(p, lt)
    lt ? 1.0 - p : p
end

function Dval(x, lp)
    lp ? log(x) : x
end

function DqIv(p, lp)
    lp ? exp(p) : p
end

function Dexp(x, lp)
    lp ? x : exp(x)
end

function Dlog(p, lp)
    lp ? p : log(p)
end

function DClog(p, lp)
    lp ? log1p(-p) : 1.0 - p
end

function Log1Exp(x)
    x > -log(2.0) ? log(-expm1(x)) : log1p(-exp(x))
end

function DLExp(x, lp)
    lp ? Log1Exp(x) : log1p(-x)
end

function DTval(x, lt, lp)
    lt ? Dval(x, lp) : DClog(x, lp)
end

function DTCval(x, lt, lp)
    lt ? DClog(x, lp) : Dval(x, lp)
end

function DTqIv(p, lt, lp)
    lp ? (lt ? exp(p) : -expm1(p)) : DLval(p, lt)
end

function DTCIv(p, lt, lp)
    lp ? (lt ? -expm1(p) : exp(p)) : DCval(p, lt)
end

function DTexp(x, lt, lp)
    Dexp(DLval(x, lt), lp)
end

function DTCexp(x, lt, lp)
    Dexp(DCval(x, lt), lp)
end

function DTlog(p, lt, lp)
    lt ? Dlog(p, lp) : DLExp(p, lp)
end

function DTClog(p, lt, lp)
    lt ? DLExp(p, lp) : Dlog(p, lp)
end

function DTLog(p, lt)
    lt ? p : Log1Exp(p)
end

function QP01check(p, lp)
    if (lp && p > 0) || (!lp && (p < 0 || p > 1))
        NaN
    end
end

function QP01boundaries(p, left, right, lt, lp)
    if lp
        if p > 0
            return NaN
        end
        if p == 0
            return lt ? right : left
        end
        if p == -Inf
            return lt ? left : right
        end
    else
        if p < 0 || p > 1
            return NaN
        end
        if p == 0
            return lt ? left : right
        end
        if p == 1
            return lt ? right : left
        end
    end
end

function Pbounds01(x, x_min, x_max, lt, lp)
    if x <= x_min
        return DT0(lt, lp)
    end
    if x >= x_max
        return DT1(lt, lp)
    end
end

function PboundsInf01(x, lt, lp)
    if isinf(x)
        if x > 0
            return DT1(lt, lp)
        end
        return DT0(lt, lp)
    end
end

function Dfexp(f, x)
    lp ? -0.5 * log(f) + x : exp(x) / √(f)
end

function nonint(x)
    abs(x) - round(x) > 1e-7 * max(1.0, abs(x))
end

function DnegInonint(x)
    x < 0.0 || nonint(x)
end

function Dnonintcheck(x, lp)
    if nonint(x)
        println("non-integer x = ", x)
        return D0(lp)
    end
end

# test session
# println(Dnonintcheck(12, false)) # return nothing
