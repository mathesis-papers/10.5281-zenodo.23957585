import argparse
import hashlib
import json
from math import factorial, isqrt
from pathlib import Path


SIGNED_PREFIX = [
    0, 1, 2, 3, 4, 1, 2, 3, 4, 5, 2, 3, 4, 5, 6, 3, 4, 5, 6, 7,
    4, 5, 6, 7, 4, 5, 6, 7, 8, 5, 6, 7, 8, 5, 6, 7, 8, 9, 6, 7,
    8, 9, 6, 7, 8, 9, 10, 7, 8, 9, 10, 7, 8, 9, 10, 11, 8, 9, 10,
    11, 8, 9, 10, 11, 12, 9, 10, 11, 12, 9, 10, 11,
]
ABSOLUTE_PREFIX = [
    0, 1, 2, 1, 2, 3, 2, 3, 2, 3, 4, 3, 4, 3, 4, 5, 4, 5, 4, 5,
    6, 5, 6, 5, 6, 5, 6, 7, 6, 7, 6, 7, 6, 7, 6, 7, 8, 7, 8, 7,
    8, 7, 8, 7, 8, 9, 8, 9, 8, 9, 8, 9, 8, 9, 10, 9, 10, 9, 10,
    9, 10, 9, 10, 9, 10, 11, 10, 11,
]
TIE_PREFIX = [
    8, 13, 34, 43, 76, 89, 134, 151, 208, 229, 298, 323, 404,
    433, 526, 559, 664, 701, 818, 859, 988, 1033, 1174, 1223,
    1376, 1429, 1594, 1651, 1828, 1889, 2078, 2143, 2344,
    2413, 2626, 2699, 2924,
]


def require(condition, label):
    if not condition:
        raise AssertionError(label)


def signed_mode(degree):
    return (2 * degree + 6 - isqrt(8 * degree + 21)) // 8


def signed_exponent(degree):
    return degree - 4 * signed_mode(degree)


def absolute_exponent(degree):
    exponent = degree % 2
    while exponent * exponent + 5 * exponent + 2 < 2 * degree:
        exponent += 2
    return exponent


def absolute_mode_closed(degree):
    discriminant = 8 * degree + 17
    root = isqrt(discriminant)
    correction = int(root * root != discriminant)
    return (2 * degree + 5 - root - correction) // 4


def tie_degree(index):
    return 2 * index * (index + 2) + (index + 1) * (-1) ** (index + 1)


def fourth(exponent):
    return (exponent + 1) * (exponent + 2) * (exponent + 3) * (exponent + 4)


def denominator(mode):
    return 16 * (2 * mode) * (2 * mode - 1)


def coefficient_checks(limit):
    older = [0]
    coefficients = [1]
    previous_absolute = 0
    found_ties = []
    for degree in range(limit + 1):
        if degree == 1:
            older, coefficients = coefficients, [0, 2]
        elif degree >= 2:
            following = [0] * (degree + 1)
            for exponent, value in enumerate(coefficients):
                following[exponent + 1] += 2 * value
            for exponent, value in enumerate(older):
                following[exponent] -= 2 * (degree - 1) * value
            older, coefficients = coefficients, following
        largest = max(coefficients)
        signed_positions = [exponent for exponent, value in enumerate(coefficients) if value == largest]
        selected = signed_exponent(degree)
        require(signed_positions == [selected], ('signed maximum', degree))
        mode = signed_mode(degree)
        value = factorial(degree) * 2 ** selected // (factorial(2 * mode) * factorial(selected))
        require(value == largest, ('signed maximum value', degree))
        for exponent, value in enumerate(coefficients):
            if (degree - exponent) % 2:
                expected = 0
            else:
                index = (degree - exponent) // 2
                expected = (-1) ** index * factorial(degree) * 2 ** exponent
                expected //= factorial(index) * factorial(exponent)
            require(value == expected, ('finite expansion', degree, exponent))
        maximum_absolute = max(map(abs, coefficients))
        absolute_positions = [exponent for exponent, value in enumerate(coefficients) if abs(value) == maximum_absolute]
        absolute = absolute_exponent(degree)
        require(absolute_positions[0] == absolute, ('absolute maximum', degree))
        require(abs(selected - absolute) <= 2, ('signed absolute interlacing', degree))
        absolute_mode = (degree - absolute) // 2
        require(absolute_mode == absolute_mode_closed(degree), ('absolute closed formula', degree))
        require((largest == maximum_absolute) == (len(absolute_positions) == 2 or absolute_mode % 2 == 0), ('signed absolute equality', degree))
        if degree:
            require(absolute - previous_absolute in (1, -1), ('absolute transition', degree))
        previous_absolute = absolute
        if len(absolute_positions) == 2:
            require(absolute_positions[1] == absolute + 2, ('adjacent tie', degree))
            found_ties.append(degree)
        else:
            require(len(absolute_positions) == 1, ('absolute multiplicity', degree))
        value = factorial(degree) * 2 ** absolute
        value //= factorial((degree - absolute) // 2) * factorial(absolute)
        require(value == maximum_absolute, ('absolute maximum value', degree))
    expected_ties = []
    index = 1
    while tie_degree(index) <= limit:
        expected_ties.append(tie_degree(index))
        index += 1
    require(found_ties == expected_ties, 'complete tie list in polynomial range')
    return len(found_ties)


def absolute_threshold_checks(limit):
    mode = 0
    previous_transition = 0
    block_count = 0
    long_blocks = 0
    for degree in range(limit + 1):
        previous_mode = mode
        candidate_exponent = degree - 2 * (mode + 1)
        if candidate_exponent >= 0 and (candidate_exponent + 1) * (candidate_exponent + 2) >= 4 * (mode + 1):
            mode += 1
        require(mode == absolute_mode_closed(degree), ('absolute threshold closed formula', degree))
        exponent = degree - 2 * mode
        require(exponent >= 0 and exponent * exponent + 5 * exponent + 2 >= 2 * degree, ('absolute threshold', degree))
        if exponent >= 2:
            lower = exponent - 2
            require(lower * lower + 5 * lower + 2 < 2 * degree, ('least absolute threshold', degree))
        if mode != previous_mode:
            length = degree - previous_transition
            require(length in (2, 3), ('absolute block length', degree))
            if previous_mode == 0:
                require(length == 3, 'first absolute block')
                expected_long = True
            else:
                threshold_exponent = previous_transition - 2 * previous_mode
                last_mode = (threshold_exponent + 1) * (threshold_exponent + 2) // 4
                expected_long = previous_mode == last_mode
                if expected_long:
                    require(previous_transition == 2 * last_mode + threshold_exponent, ('absolute long start', degree))
            require((length == 3) == expected_long, ('absolute long block indices', degree))
            block_count += 1
            long_blocks += int(length == 3)
            previous_transition = degree
    return block_count, long_blocks


def value_ratio_checks(limit):
    previous_signed = 1
    previous_absolute = 1
    for degree in range(1, limit + 1):
        earlier = degree - 1
        signed_exponent_before = signed_exponent(earlier)
        signed_mode_before = signed_mode(earlier)
        selected = signed_exponent(degree)
        signed_value = factorial(degree) * 2 ** selected // (factorial(2 * signed_mode(degree)) * factorial(selected))
        if signed_mode(degree) == signed_mode_before:
            numerator = 2 * degree
            divisor = signed_exponent_before + 1
        else:
            numerator = degree * signed_exponent_before * (signed_exponent_before - 1) * (signed_exponent_before - 2)
            divisor = 8 * (2 * signed_mode_before + 2) * (2 * signed_mode_before + 1)
        require(signed_value * divisor == previous_signed * numerator, ('signed exact ratio', degree))
        absolute_before = absolute_exponent(earlier)
        mode_before = (earlier - absolute_before) // 2
        absolute = absolute_exponent(degree)
        selected_mode = (degree - absolute) // 2
        absolute_value = factorial(degree) * 2 ** absolute // (factorial(selected_mode) * factorial(absolute))
        if selected_mode == mode_before:
            numerator = 2 * degree
            divisor = absolute_before + 1
        else:
            numerator = degree * absolute_before
            divisor = 2 * (mode_before + 1)
        require(absolute_value * divisor == previous_absolute * numerator, ('absolute exact ratio', degree))
        previous_signed, previous_absolute = signed_value, absolute_value
    return limit


def threshold_checks(limit):
    previous_mode = 0
    threshold_mode = 0
    previous_transition = 0
    transition_count = 0
    long_blocks = 0
    for degree in range(limit + 1):
        while 4 * (threshold_mode + 1) <= degree:
            exponent = degree - 4 * (threshold_mode + 1)
            if fourth(exponent) < denominator(threshold_mode + 1):
                break
            threshold_mode += 1
        mode = signed_mode(degree)
        require(mode == threshold_mode, ('quartic threshold', degree))
        exponent = degree - 4 * mode
        require(0 <= exponent <= degree, ('support', degree))
        require(exponent * exponent + 7 * exponent + 7 >= 2 * degree, ('quadratic threshold', degree))
        if exponent >= 4:
            lower = exponent - 4
            require(lower * lower + 7 * lower + 7 < 2 * degree, ('minimal quadratic threshold', degree))
        if degree:
            change = mode - previous_mode
            require(change in (0, 1), ('mode transition', degree))
            if change:
                block_length = degree - previous_transition
                require(block_length in (4, 5), ('block length', degree))
                old_mode = mode - 1
                if old_mode == 0:
                    require(block_length == 5, 'initial block')
                    expected_long = True
                else:
                    threshold_exponent = previous_transition - 4 * old_mode
                    last_mode = (threshold_exponent ** 2 + 5 * threshold_exponent + 7) // 8
                    expected_long = old_mode == last_mode
                    require(threshold_exponent ** 2 + 5 * threshold_exponent + 7 >= 8 * old_mode, ('block start threshold', degree))
                require((block_length == 5) == expected_long, ('quadratic long-block indices', degree))
                transition_count += 1
                long_blocks += int(block_length == 5)
                previous_transition = degree
        previous_mode = mode
    return transition_count, long_blocks


def generating_function_checks(limit):
    values = [0] + [tie_degree(index) for index in range(1, limit + 1)]
    numerator = [0, 8, 5, 5, -1, -1]
    divisor = [1, -1, -2, 2, 1, -1]
    for index in range(limit + 1):
        product = sum(divisor[offset] * values[index - offset] for offset in range(min(index, 5) + 1))
        require(product == (numerator[index] if index < len(numerator) else 0), ('generating function', index))
    return limit


def main():
    parser = argparse.ArgumentParser(description='Exact finite checks for Hermite coefficient maxima.')
    parser.add_argument('--output', type=Path)
    arguments = parser.parse_args()
    polynomial_limit = 600
    threshold_limit = 1_000_000
    generating_limit = 10_000
    require([signed_exponent(degree) for degree in range(len(SIGNED_PREFIX))] == SIGNED_PREFIX, 'A381423 displayed prefix')
    require([absolute_exponent(degree) for degree in range(len(ABSOLUTE_PREFIX))] == ABSOLUTE_PREFIX, 'A381524 displayed prefix')
    require([tie_degree(index) for index in range(1, len(TIE_PREFIX) + 1)] == TIE_PREFIX, 'A381737 displayed prefix')
    ties = coefficient_checks(polynomial_limit)
    transitions, long_blocks = threshold_checks(threshold_limit)
    absolute_blocks, absolute_long_blocks = absolute_threshold_checks(threshold_limit)
    value_ratio_checks(polynomial_limit)
    generating_function_checks(generating_limit)
    result = {
        'status': 'pass',
        'arithmetic': 'exact integers',
        'coefficient_recurrence': 'H_0=1; H_1=2x; H_n=2xH_(n-1)-2(n-1)H_(n-2)',
        'polynomial_degrees': {'first': 0, 'last': polynomial_limit, 'count': polynomial_limit + 1},
        'polynomial_checks': ['finite expansion', 'unique signed maximum', 'signed maximum value', 'least absolute maximum', 'absolute closed formula', 'absolute difference values', 'absolute maximum value', 'complete two-maximum degree list', 'signed absolute equality criterion', 'signed absolute interlacing'],
        'two_maximum_degrees_in_polynomial_range': ties,
        'threshold_degrees': {'first': 0, 'last': threshold_limit, 'count': threshold_limit + 1},
        'threshold_checks': ['quartic and quadratic selectors', 'closed formula', 'mode increments', 'difference-block lengths', 'quadratic indices of length-five blocks'],
        'complete_blocks': transitions,
        'length_five_blocks': long_blocks,
        'absolute_threshold_checks': ['closed formula with exact-square boundary', 'least quadratic selector', 'two difference-block lengths', 'quadratic length-three block indices'],
        'absolute_complete_blocks': absolute_blocks,
        'absolute_length_three_blocks': absolute_long_blocks,
        'exact_value_ratio_degrees': {'first': 1, 'last': polynomial_limit, 'count': polynomial_limit},
        'generating_function_coefficients': {'first': 0, 'last': generating_limit, 'count': generating_limit + 1},
        'oeis_prefixes': {'A381423': len(SIGNED_PREFIX), 'A381524': len(ABSOLUTE_PREFIX), 'A381737': len(TIE_PREFIX)},
        'scope': 'Finite comparisons; the all-degree conclusions have written proofs in note.tex.',
        'checker_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    text = json.dumps(result, indent=2) + '\n'
    if arguments.output:
        arguments.output.write_text(text)
    print(text, end='')


if __name__ == '__main__':
    main()
