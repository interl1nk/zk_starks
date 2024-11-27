from hashlib import sha256
from random import randint

from fields.field import FieldElement
from merkle import MerkleTree, verify_decommitment


def test_merkle_get_authentication_path():
    for _ in range(10):
        data_length = randint(0, 2000)
        data = [FieldElement.random_element() for _ in range(data_length)]
        m = MerkleTree(data)
        leaf_id = randint(0, data_length - 1)
        decommitment = m.get_authentication_path(leaf_id)
        # Проверьте правильность декоммита.
        content = data[leaf_id]
        assert verify_decommitment(leaf_id, content, decommitment, m.root)
        # Проверьте, что изменение декоммитента приводит к сбою проверки.
        altered = decommitment[:]
        random_index = randint(0, len(altered) - 1)
        altered[random_index] = sha256(altered[random_index].encode()).hexdigest()
        assert not verify_decommitment(leaf_id, content, altered, m.root)
        # Проверьте, что изменение содержимого приводит к сбою проверки.
        other_content = data[randint(0, data_length - 1)]
        assert not verify_decommitment(
            leaf_id, other_content, decommitment, m.root) or other_content == content