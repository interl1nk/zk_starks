import operator
try:
    from tqdm import tqdm
except ModuleNotFoundError:
    # `tqdm` - это обертка для итераторов, реализующая прогресс-бар. Если это
    # недоступно, просто верните сам итератор.
    tqdm = lambda x: x

from fields.field import FieldElement
from utils.list_utils import remove_trailing_elements, scalar_operation, two_lists_tuple_operation


def trim_trailing_zeros(p):
    """
    Удаляет нули из конца списка.
    """
    return remove_trailing_elements(p, FieldElement.zero())


def prod(values):
    """
    Вычисляет продукт.
    """
    len_values = len(values)
    if len_values == 0:
        return 1
    if len_values == 1:
        return values[0]
    return prod(values[:len_values // 2]) * prod(values[len_values // 2:])


def latex_monomial(exponent, coef, var):
    """
    Возвращает строковое представление монома в формате LaTeX.
    """
    if exponent == 0:
        return str(coef)
    if coef == 1:
        coef = ''
    if coef == -1:
        coef = '-'
    if exponent == 1:
        return f'{coef}{var}'
    return f'{coef}{var}^{{{exponent}}}'


class Polynomial:
    """
    Представляет полином над FieldElement.
    """

    @classmethod
    def X(cls):
        """
        Возвращает многочлен x.
        """
        return cls([FieldElement.zero(), FieldElement.one()])

    def __init__(self, coefficients, var='x'):
        # Внутреннее хранение коэффициентов в self.poly, наименее значимый (т.е. свободный член)
        # первым, поэтому $9 - 3x^2 + 19x^5$ представляется внутри списка [9, 0, -3, 0, 0, 19].
        # Обратите внимание, что коэффициенты копируются, поэтому вызывающая сторона может свободно изменять данный аргумент.
        self.poly = remove_trailing_elements(coefficients, FieldElement.zero())
        self.var = var

    def _repr_latex_(self):
        """
        Возвращает представление полинома в формате LaTeX для Jupyter.
        """
        if not self.poly:
            return '$0$'
        res = ['$']
        first = True
        for exponent, coef in enumerate(self.poly):
            if coef == 0:
                continue
            monomial = latex_monomial(exponent, coef, self.var)
            if first:
                first = False
                res.append(monomial)
                continue
            oper = '+'
            if monomial[0] == '-':
                oper = '-'
                monomial = monomial[1:]
            res.append(oper)
            res.append(monomial)
        res.append('$')
        return ' '.join(res)

    def __eq__(self, other):
        try:
            other = Polynomial.typecast(other)
        except AssertionError:
            return False
        return self.poly == other.poly

    @staticmethod
    def typecast(other):
        """
        Создает многочлен из `FieldElement` или `int`.
        """
        if isinstance(other, int):
            other = FieldElement(other)
        if isinstance(other, FieldElement):
            other = Polynomial([other])
        assert isinstance(other, Polynomial), f'Type mismatch: Polynomial and {type(other)}.'
        return other

    def __add__(self, other):
        other = Polynomial.typecast(other)
        return Polynomial(two_lists_tuple_operation(
            self.poly, other.poly, operator.add, FieldElement.zero()))

    __radd__ = __add__  # Для поддержки <int> + <Polynomial> (как в `1 + x + x**2`).

    def __sub__(self, other):
        other = Polynomial.typecast(other)
        return Polynomial(two_lists_tuple_operation(
            self.poly, other.poly, operator.sub, FieldElement.zero()))

    def __rsub__(self, other):  # Для поддержки <int> - <Polynomial> (как в `1 - x + x**2`).
        return -(self - other)

    def __neg__(self):
        return Polynomial([]) - self

    def __mul__(self, other):
        other = Polynomial.typecast(other)
        pol1, pol2 = [[x.val for x in p.poly] for p in (self, other)]
        res = [0] * (self.degree() + other.degree() + 1)
        for i, c1 in enumerate(pol1):
            for j, c2 in enumerate(pol2):
                res[i + j] += c1 * c2
        res = [FieldElement(x) for x in res]
        return Polynomial(res)

    __rmul__ = __mul__  # Для поддержки <int> * <Polynomial>.

    def compose(self, other):
        """
        Составляет данный многочлен с `другими`.
        Пример:
        >>> f = X**2 + X
        >>> g = X + 1
        >>> f.compose(g) == (2 + 3*X + X**2)
        True
        """
        other = Polynomial.typecast(other)
        res = Polynomial([])
        for coef in self.poly[::-1]:
            res = (res * other) + Polynomial([coef])
        return res

    def qdiv(self, other):
        """
        Возвращает q, r соответственно многочлены с коэффициентом и остатком, такие что
        f = q * g + r, где deg(r) < deg(g).
        * Утверждение, что g не является нулевым многочленом.
        """
        other = Polynomial.typecast(other)
        pol2 = trim_trailing_zeros(other.poly)
        assert pol2, 'Деление на нулевой многочлен.'
        pol1 = trim_trailing_zeros(self.poly)
        if not pol1:
            return [], []
        rem = pol1
        deg_dif = len(rem) - len(pol2)
        quotient = [FieldElement.zero()] * (deg_dif + 1)
        g_msc_inv = pol2[-1].inverse()
        while deg_dif >= 0:
            tmp = rem[-1] * g_msc_inv
            quotient[deg_dif] = quotient[deg_dif] + tmp
            last_non_zero = deg_dif - 1
            for i, coef in enumerate(pol2, deg_dif):
                rem[i] = rem[i] - (tmp * coef)
                if rem[i] != FieldElement.zero():
                    last_non_zero = i
            # Устраните нули в конце (т.е. сделайте так, чтобы r заканчивался последним ненулевым коэффициентом).
            rem = rem[:last_non_zero + 1]
            deg_dif = len(rem) - len(pol2)
        return Polynomial(trim_trailing_zeros(quotient)), Polynomial(rem)

    def __truediv__(self, other):
        div, mod = self.qdiv(other)
        assert mod == 0, 'Многочлены не являются делимыми.'
        return div

    def __mod__(self, other):
        return self.qdiv(other)[1]

    @staticmethod
    def monomial(degree, coefficient):
        """
        Строит коэффициент монома * x**degree.
        """
        return Polynomial([FieldElement.zero()] * degree + [coefficient])

    @staticmethod
    def gen_linear_term(point):
        """
        Генерирует многочлен (x-p) для заданной точки p.
        """
        return Polynomial([FieldElement.zero() - point, FieldElement.one()])

    def degree(self):
        """
        Многочлены представлены в виде списка, поэтому степень равна длине списка минус
        количество нулей в конце списка (если они есть) минус 1.
        Это означает, что степень нулевого многочлена будет равна -1.
        """
        return len(trim_trailing_zeros(self.poly)) - 1

    def get_nth_degree_coefficient(self, n):
        """
        Возвращает коэффициент x**n
        """
        if n > self.degree():
            return FieldElement.zero()
        else:
            return self.poly[n]

    def scalar_mul(self, scalar):
        """
        Умножает многочлен на скаляр
        """
        return Polynomial(scalar_operation(self.poly, operator.mul, scalar))

    def eval(self, point):
        """
        Оценивает многочлен в заданной точке с помощью оценки Хорнера.
        """
        point = FieldElement.typecast(point).val
        # Работа с цифрами (в отличие от `FieldElement`) значительно ускоряет выполнение eval.
        val = 0
        for coef in self.poly[::-1]:
            val = (val * point + coef.val) % FieldElement.k_modulus
        return FieldElement(val)

    def __call__(self, other):
        """
        Если `other` - int или FieldElement, вычисляется многочлен на `other` (в поле).
        Если `other` является многочленом, то компонует self с `other` как self(other(x)).
        """
        if isinstance(other, (int)):
            other = FieldElement(other)
        if isinstance(other, FieldElement):
            return self.eval(other)
        if isinstance(other, Polynomial):
            return self.compose(other)
        raise NotImplementedError()

    def __pow__(self, other):
        """
        Вычисляет self**other с помощью повторного возведения в квадрат.
        """
        assert other >= 0
        res = Polynomial([FieldElement(1)])
        cur = self
        while True:
            if other % 2 != 0:
                res *= cur
            other >>= 1
            if other == 0:
                break
            cur = cur * cur
        return res


# python представление формальной переменной x.
X = Polynomial.X()


def calculate_lagrange_polynomials(x_values):
    """
    Учитывая значения x_ для оценки некоторых полиномов, он вычисляет часть полиномов Лагранжа
    необходимых для интерполяции полинома над этой областью.
    """
    lagrange_polynomials = []
    monomials = [Polynomial.monomial(1, FieldElement.one()) -
                 Polynomial.monomial(0, x) for x in x_values]
    numerator = prod(monomials)
    for j in tqdm(range(len(x_values))):
        # В знаменателе имеем:
        # (x_j-x_0)(x_j-x_1)...(x_j-x_{j-1})(x_j-x_{j+1})...(x_j-x_{len(X)-1})
        denominator = prod([x_values[j] - x for i, x in enumerate(x_values) if i != j])
        # Числитель немного сложнее, так как здесь нам нужно вычислить полиумножение.
        # Аналогично со знаменателем, имеем:
        # (x-x_0)(x-x_1)...(x-x_{j-1})(x-x_{j+1})...(x-x_{len(X)-1})
        cur_poly, _ = numerator.qdiv(monomials[j].scalar_mul(denominator))
        lagrange_polynomials.append(cur_poly)
    return lagrange_polynomials


def interpolate_poly_lagrange(y_values, lagrange_polynomials):
    """
    :param y_values: y координаты y точек.
    :param lagrange_polynomials: полиномы, полученные в результате calculate_lagrange_polynomials.
    :return: 'the interpolated poly'
    """
    poly = Polynomial([])
    for j, y_value in enumerate(y_values):
        poly += lagrange_polynomials[j].scalar_mul(y_value)
    return poly


def interpolate_poly(x_values, y_values):
    """
    Возвращает многочлен степени < len(x_values), который оценивает y_values[i] по x_values[i] для
    всех i.
    """
    assert len(x_values) == len(y_values)
    assert all(isinstance(val, FieldElement) for val in x_values),\
        'Не все значения x_values являются FieldElement'
    lp = calculate_lagrange_polynomials(x_values)
    assert all(isinstance(val, FieldElement) for val in y_values),\
        'Не все значения y_values являются FieldElement'
    return interpolate_poly_lagrange(y_values, lp)
