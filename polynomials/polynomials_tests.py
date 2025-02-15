from random import randint

from fields.field import FieldElement
from polynomial import Polynomial, interpolate_poly, X, prod


def random_polynomial(degree):
    """
    Возвращает случайный многочлен заданной степени, который НЕ является нулевым многочленом.
    """
    leading = FieldElement.random_element(exclude_elements=[FieldElement.zero()])
    p = [FieldElement.random_element() for _ in range(degree)] + [leading]
    return Polynomial(p)


def test_div_rand_poly():
    """
    Пусть a, b - случайные многочлены.
    Пусть q, r - это купон и остаток, которые получаются в результате выполнения функции
    qdiv(a, b).
    Мы проверяем, что действительно a = q * b + r.
    """
    iterations = 20
    for _ in range(iterations):
        deg_a = randint(0, 50)
        deg_b = randint(0, 50)
        a = random_polynomial(deg_a)
        b = random_polynomial(deg_b)
        (q, r) = a.qdiv(b)
        d = r + q * b
        assert r.degree() < b.degree()
        assert d == a


def test_poly_interpolation():
    for _ in range(10):
        # Сгенерируйте случайный многочлен.
        degree = randint(0, 100)
        p = random_polynomial(degree)
        # Оцените его по числу точек, которое не меньше его степени.
        x_vals = set()
        while len(x_vals) < degree + 1:
            x_vals.add(FieldElement.random_element())
        x_vals = list(x_vals)
        y_vals = [p.eval(x) for x in x_vals]
        # Получите многочлен из оценки.
        interpolated_p = interpolate_poly(x_vals, y_vals)
        # Проверьте равенство.
        assert p == interpolated_p


def test_compose():
    for _ in range(10):
        outer_poly = random_polynomial(randint(0, 1024))
        inner_poly = random_polynomial(randint(0, 16))
        # Проверьте оценку композиции poly outer_poly(inner_poly) на случайной точке.
        point_to_eval = FieldElement.random_element()
        assert((outer_poly.compose(inner_poly)).eval(point_to_eval) ==
               outer_poly.eval(inner_poly.eval(point_to_eval)))


def test_latex_repr():
    assert (X - 5 * X**2)._repr_latex_() == '$ x - 5x^{2} $'
    assert (-71 * X**5 + X**2048)._repr_latex_() == '$ -71x^{5} + x^{2048} $'
    assert (X**2 + X + 1)._repr_latex_() == '$ 1 + x + x^{2} $'
    assert (92 * X**65 + 4 * X**15 - 31)._repr_latex_() == '$ -31 + 4x^{15} + 92x^{65} $'
    assert (-X)._repr_latex_() == '$ -x $'
    assert (X-X)._repr_latex_() == '$0$'


def test_poly_mul():
    assert (X+1) * (X+1) == X**2 + 2*X + 1


def test_prod():
    g = FieldElement.generator()**((FieldElement.k_modulus - 1) // 1024)
    assert X**1024 - 1 == prod([(X - g**i) for i in range(1024)])


def test_call_compose():
    p = X**2 + X
    assert p(X + 1) == X**2 + 3*X + 2


def test_call_eval():
    p = X**2 + X
    assert p(5) == 30


def test_truediv():
    p = X**2 - 1
    assert p / (X - 1) == X + 1


def test_mod():
    p = X**9 - 5*X + 4
    assert p % (X**2 + 1) == -4*X + 4
