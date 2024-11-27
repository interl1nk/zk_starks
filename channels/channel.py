import inspect
from hashlib import sha256

from fields.field import FieldElement


def serialize(obj):
    """
    Сериализует объект в строку.
    """
    if isinstance(obj, (list, tuple)):
        return ','.join(map(serialize, obj))
    return obj._serialize_()


class Channel(object):
    """
    Экземпляр канала может быть использован проверяющим или верификатором для сохранения семантики
    интерактивной системы доказательств, хотя на самом деле он неинтерактивен и использует Sha256
    для генерации случайности, когда это требуется.
    Он позволяет записывать в него данные в строковой форме и считывать из него случайные целые числа или случайные
    FieldElements.
    """

    def __init__(self):
        self.state = '0'
        self.proof = []

    def send(self, s):
        self.state = sha256((self.state + s).encode()).hexdigest()
        self.proof.append(f'{inspect.stack()[0][3]}:{s}')

    def receive_random_int(self, min, max, show_in_proof=True):
        """
        Эмулирует случайное целое число, отправленное верификатором в диапазоне [min, max] (включая min и
        max).
        """

        # Обратите внимание, что когда диапазон близок к 2^256, это не приводит к равномерному распределению,
        # даже если sha256 распределен равномерно.
        # Однако для целей данного руководства это достаточно близко.
        num = min + (int(self.state, 16) % (max - min + 1))
        self.state = sha256((self.state).encode()).hexdigest()
        if show_in_proof:
            self.proof.append(f'{inspect.stack()[0][3]}:{num}')
        return num

    def receive_random_field_element(self):
        """
        Эмулирует элемент случайного поля, отправленный верификатором.
        """
        num = self.receive_random_int(0, FieldElement.k_modulus - 1, show_in_proof=False)
        self.proof.append(f'{inspect.stack()[0][3]}:{num}')
        return FieldElement(num)
