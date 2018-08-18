from numpy import sqrt, arctan2, linspace, meshgrid


def σ(x1, y1, x2, y2, P, a, λ, x, y):
    c = (x2 - x1) / sqrt((x2 - x1)**2 + (y2 - y1)**2)
    s = (y2 - y1) / sqrt((x2 - x1)**2 + (y2 - y1)**2)
    term0 = (x*c + y*s)
    term1 = (x1*c - y1*s) - term0
    term2 = (-x*s + y*c)
    term3 = (-x2*s + y2*c) - term2
    term4 = (x2*c + y2*s) - term0

    c = -P * a**2 * λ
    xx_ = c * (term1 / (term1 + term3)**2 - term4 / (term4 + term3)**2)
    xy_ = c * (term3 / (term3 + term1)**2 - term3 / (term3 + term4)**2)

    xx = (c**2 - s**2) * xx_ - 2*c*s*xy_
    xy = 2*c*s*xx_ + (c**2 - s**2) * xy_
    return xx, xy


def ret(σ_xx, σ_xy, C, L):
    return C * L * sqrt((2 * σ_xx)**2 + 4*σ_xy**2)


def Θ(σ_xx, σ_xy):
    return arctan2(2 * σ_xy, 2 * σ_xx)


if __name__ == '__main__':
    x, y = meshgrid(linspace(-2000, 2000, num=20),
                    linspace(-4000, 4000, num=40))
    σ_xx, σ_yy = σ(x1=0, y1=-1000, x2=0, y2=1000, P=5, a=1, λ=1/10, x=x, y=y)
    ret(σ_xx, σ_yy, C=3, L=1)
    Θ(σ_xx, σ_yy)
