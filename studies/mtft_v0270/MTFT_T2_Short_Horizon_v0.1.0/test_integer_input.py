"""Regression control for lossy JSON integer intake; standard library only."""
import copy
from dyadic_bounds import decode


def main():
    base={'exponent':0,'real':[['1']],'imag':[['0']]}
    for value in [1.5, 1.0, True, '1.5']:
        bad=copy.deepcopy(base);bad['real'][0][0]=value
        try:
            decode(bad)
        except ValueError:
            continue
        raise AssertionError(f'Lossy integer intake accepted {value!r}')
    assert decode(base)==([[1]],[[0]])
    print('PASS: four lossy integer inputs rejected; exact integer input accepted')


if __name__=='__main__':
    main()
