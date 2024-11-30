from field import FieldElement


def test_field_operations():
    # Проверьте pow, mul и модульные операции.
    t = FieldElement(2).pow(30) * FieldElement(3) + FieldElement(1)
    assert t == FieldElement(0)
    # Проверка генератора
    # Проверка инверсии
    # Проверка хэша через использование множества


def test_field_div():
    for _ in range(100):
        t = FieldElement.random_element(exclude_elements=[FieldElement.zero()])
        t_inv = FieldElement.one() / t
        assert t_inv == t.inverse()
        assert t_inv * t == FieldElement.one()
