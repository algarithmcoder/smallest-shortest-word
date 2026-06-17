import re
from typing import List, Dict, Set, Tuple
from functools import lru_cache

# ===================== Core Constants =====================
INVERSE = {'a': 'A', 'A': 'a', 'b': 'B', 'B': 'b'}
PERMUTATIONS = {
    'P1': {'a': 'a', 'b': 'b'},
    'P2': {'a': 'a', 'b': 'B'},
    'P3': {'a': 'A', 'b': 'b'},
    'P4': {'a': 'A', 'b': 'B'},
    'P5': {'a': 'b', 'b': 'a'},
    'P6': {'a': 'b', 'b': 'A'},
    'P7': {'a': 'B', 'b': 'a'},
    'P8': {'a': 'B', 'b': 'A'}
}
DEHN_TWISTS = {
    'T1': {'a': 'a', 'b': ['a', 'b']},
    'T2': {'a': 'a', 'b': ['A', 'b']},
    'T3': {'a': ['b', 'a'], 'b': 'b'},
    'T4': {'a': ['B', 'a'], 'b': 'b'}
}


# ===================== Core Functions =====================
@lru_cache(maxsize=None)
def parse_group_element(element_str):
    s = element_str.replace(' ', '').strip()
    original_s = s

    s = re.sub(r'A(?!\^)', 'a^-1', s)
    s = re.sub(r'B(?!\^)', 'b^-1', s)
    s = re.sub(r'A\^([+-]?\d+)', r'a^-1^\1', s)
    s = re.sub(r'B\^([+-]?\d+)', r'b^-1^\1', s)

    if not s:
        raise ValueError("Input cannot be empty")

    parsed = []
    i = 0
    n = len(s)

    def parse_exponent(start_idx):
        idx = start_idx
        exp_str = ""
        if idx < n and s[idx] in '+-':
            exp_str += s[idx]
            idx += 1
        while idx < n and s[idx].isdigit():
            exp_str += s[idx]
            idx += 1
        if not exp_str or (len(exp_str) == 1 and exp_str[0] in '+-'):
            raise ValueError(f"Invalid exponent: {exp_str}")
        return int(exp_str), idx

    while i < n:
        if s[i] == ' ':
            i += 1
            continue
        if s[i] == '(':
            depth = 1
            j = i + 1
            while j < n and depth > 0:
                if s[j] == '(':
                    depth += 1
                elif s[j] == ')':
                    depth -= 1
                j += 1
            if depth > 0:
                raise ValueError(f"Mismatched parentheses: {original_s}")
            content = s[i + 1:j - 1]
            i = j - 1
            if j < n and s[j] == '^':
                j += 1
                exp, j = parse_exponent(j)
                i = j
            else:
                exp = 1
                i = j
            if content:
                sub_parsed = parse_group_element(content)
                if exp == 0:
                    pass
                elif exp > 0:
                    parsed.extend(sub_parsed * exp)
                else:
                    reversed_sub = [(g, -p) for g, p in reversed(sub_parsed)]
                    parsed.extend(reversed_sub * (-exp))
        elif s[i] in 'ab':
            gen = s[i]
            i += 1
            total_exp = 1
            while i < n and s[i] == '^':
                i += 1
                exp, i = parse_exponent(i)
                total_exp *= exp
            if total_exp != 0:
                parsed.append((gen, total_exp))
        else:
            if i + 2 < n and s[i:i + 3] in ['a^-', 'b^-']:
                gen = s[i]
                i += 3
                exp, i = parse_exponent(i)
                exp = -exp
                while i < n and s[i] == '^':
                    i += 1
                    next_exp, i = parse_exponent(i)
                    exp *= next_exp
                if exp != 0:
                    parsed.append((gen, exp))
            else:
                raise ValueError(f"Unrecognized character: '{s[i]}' at position {i}, only a, b, A, B, parentheses() and exponents^n are supported")
    return parsed


def reduce_redundancy(parsed_units):
    if not parsed_units:
        return []
    reduced = []
    current_gen, current_power = parsed_units[0]
    for gen, power in parsed_units[1:]:
        if gen == current_gen:
            current_power += power
        else:
            if current_power != 0:
                reduced.append((current_gen, current_power))
            current_gen, current_power = gen, power
    if current_power != 0:
        reduced.append((current_gen, current_power))
    return reduced


def expand_to_letters(reduced_units):
    letters = []
    for gen, power in reduced_units:
        if power > 0:
            letters.extend([gen] * power)
        else:
            inv_gen = 'A' if gen == 'a' else 'B'
            letters.extend([inv_gen] * (-power))
    return letters


def reduce_adjacent_letters(letters):
    stack = []
    for letter in letters:
        if stack:
            top = stack[-1]
            if (top == 'a' and letter == 'A') or (top == 'A' and letter == 'a'):
                stack.pop()
            elif (top == 'b' and letter == 'B') or (top == 'B' and letter == 'b'):
                stack.pop()
            else:
                stack.append(letter)
        else:
            stack.append(letter)
    return stack


def letters_to_reduced(letters):
    clean_letters = reduce_adjacent_letters(letters)
    if not clean_letters:
        return []
    parsed = []
    for letter in clean_letters:
        if letter == 'a':
            parsed.append(('a', 1))
        elif letter == 'A':
            parsed.append(('a', -1))
        elif letter == 'b':
            parsed.append(('b', 1))
        elif letter == 'B':
            parsed.append(('b', -1))
    return reduce_redundancy(parsed)


def get_cyclic_shifts(letters):
    shifts = []
    n = len(letters)
    for i in range(n):
        shifted = letters[i:] + letters[:i]
        shifts.append(shifted)
    return shifts


def calculate_cyclic_reduced_length(reduced_units):
    if not reduced_units:
        return 0
    letters = expand_to_letters(reduced_units)
    if not letters:
        return 0
    shifts = get_cyclic_shifts(letters)
    min_length = float('inf')
    for shifted in shifts:
        reduced_shifted = letters_to_reduced(shifted)
        length = sum(abs(p) for (g, p) in reduced_shifted)
        if length < min_length:
            min_length = length
    return min_length


def apply_automorphism_1(letters):
    new_letters = []
    for letter in letters:
        if letter == 'b':
            new_letters.extend(['a', 'b'])
        elif letter == 'B':
            new_letters.extend(['B', 'A'])
        else:
            new_letters.append(letter)
    return reduce_adjacent_letters(new_letters)


def apply_automorphism_2(letters):
    new_letters = []
    for letter in letters:
        if letter == 'b':
            new_letters.extend(['A', 'b'])
        elif letter == 'B':
            new_letters.extend(['B', 'a'])
        else:
            new_letters.append(letter)
    return reduce_adjacent_letters(new_letters)


def apply_automorphism_3(letters):
    new_letters = []
    for letter in letters:
        if letter == 'a':
            new_letters.extend(['b', 'a'])
        elif letter == 'A':
            new_letters.extend(['A', 'B'])
        else:
            new_letters.append(letter)
    return reduce_adjacent_letters(new_letters)


def apply_automorphism_4(letters):
    new_letters = []
    for letter in letters:
        if letter == 'a':
            new_letters.extend(['B', 'a'])
        elif letter == 'A':
            new_letters.extend(['A', 'b'])
        else:
            new_letters.append(letter)
    return reduce_adjacent_letters(new_letters)


def format_result(reduced_units):
    if not reduced_units:
        return 'e'
    parts = []
    for gen, power in reduced_units:
        if power == 1:
            parts.append(gen)
        elif power == -1:
            parts.append('A' if gen == 'a' else 'B')
        elif power > 1:
            parts.append(f"{gen}^{power}")
        else:
            inv_gen = 'A' if gen == 'a' else 'B'
            parts.append(f"{inv_gen}^{abs(power)}")
    return ''.join(parts)


def get_shortest_cyclic_reduced_word(input_str):
    raw_units = parse_group_element(input_str)
    current_reduced = reduce_redundancy(raw_units)
    current_length = calculate_cyclic_reduced_length(current_reduced)
    current_str = format_result(current_reduced)
    processed = {current_str}
    has_improved = True

    while has_improved:
        has_improved = False
        current_letters = expand_to_letters(current_reduced)
        candidates = []
        auto_functions = [
            (apply_automorphism_1, 1),
            (apply_automorphism_2, 2),
            (apply_automorphism_3, 3),
            (apply_automorphism_4, 4)
        ]

        for auto_func, auto_id in auto_functions:
            auto_letters = auto_func(current_letters)
            auto_reduced = letters_to_reduced(auto_letters)
            auto_length = calculate_cyclic_reduced_length(auto_reduced)
            auto_str = format_result(auto_reduced)
            candidates.append((auto_length, auto_reduced, auto_str, auto_id))

        better_candidates = [c for c in candidates if c[0] < current_length]
        if better_candidates:
            better_candidates.sort(key=lambda x: (x[0], x[3]))
            best = better_candidates[0]
            best_length, best_reduced, best_str, best_auto_id = best
            if best_str not in processed:
                current_reduced = best_reduced
                current_length = best_length
                current_str = best_str
                processed.add(current_str)
                has_improved = True
    return current_str


def invert_word(chars: List[str]) -> List[str]:
    return [INVERSE[c] for c in reversed(chars)]


def normalize_input(expr: str) -> str:
    expr = expr.replace(' ', '')
    expr = re.sub(r'a\^-1', 'A', expr)
    expr = re.sub(r'b\^-1', 'B', expr)
    return expr


def expand_exponent(expr: str) -> str:
    def replace_neg(match):
        base = match.group(1)
        n = int(match.group(2))
        return INVERSE[base] * n

    def replace_pos(match):
        base = match.group(1)
        n = int(match.group(2))
        return base * n

    changed = True
    while changed:
        changed = False
        new_expr = re.sub(r'([ab])\^-(\d+)', replace_neg, expr)
        if new_expr != expr:
            expr = new_expr
            changed = True
        new_expr = re.sub(r'([abAB])\^(\d+)', replace_pos, expr)
        if new_expr != expr:
            expr = new_expr
            changed = True
    return expr


def parse_expression(expr: str) -> List[str]:
    expr = normalize_input(expr)
    while '(' in expr:
        match = re.search(r'\(([^()]+)\)\^(-?\d+)', expr)
        if match:
            inner = match.group(1)
            exp = int(match.group(2))
            expanded_inner = expand_exponent(inner)
            inner_chars = list(expanded_inner)
            if exp > 0:
                replacement = expanded_inner * exp
            elif exp < 0:
                inv_chars = invert_word(inner_chars)
                replacement = ''.join(inv_chars) * (-exp)
            else:
                replacement = ''
            expr = expr.replace(match.group(0), replacement)
        else:
            expr = re.sub(r'\(([^()]+)\)', r'\1', expr)
    expr = expand_exponent(expr)
    chars = []
    i = 0
    while i < len(expr):
        if expr[i] in 'abAB':
            chars.append(expr[i])
        i += 1
    return chars


def reduce_word(chars: List[str]) -> List[str]:
    result = []
    for c in chars:
        if result and result[-1] == INVERSE[c]:
            result.pop()
        else:
            result.append(c)
    return result


def cyclic_reduce(chars: List[str]) -> List[str]:
    if len(chars) <= 1:
        return chars.copy()
    result = chars.copy()
    while len(result) > 1 and result[0] == INVERSE[result[-1]]:
        result.pop(0)
        result.pop()
    return result


class FreeGroupWord:
    def __init__(self, chars: List[str]):
        self.raw_chars = chars
        self.reduced_chars = reduce_word(chars)
        self.cyclic_chars = cyclic_reduce(self.reduced_chars)
        self.length = len(self.cyclic_chars)
        self._validate_cyclic_reduced()

    def _validate_cyclic_reduced(self):
        temp = self.cyclic_chars.copy()
        if len(temp) <= 1:
            return
        if temp[0] == INVERSE[temp[-1]]:
            raise ValueError(f"Word {''.join(temp)} is not cyclically reduced! First and last characters {temp[0]} and {temp[-1]} are inverses")
        for i in range(len(temp) - 1):
            if temp[i] == INVERSE[temp[i + 1]]:
                raise ValueError(
                    f"Word {''.join(temp)} is not reduced! Characters at positions {i} and {i + 1} ({temp[i]} and {temp[i + 1]}) are inverses")

    def get_length(self) -> int:
        return self.length

    def __str__(self) -> str:
        return ''.join(self.cyclic_chars)

    def get_cyclic_permutations(self) -> Set[str]:
        n = len(self.cyclic_chars)
        permutations = set()
        for i in range(n):
            perm = self.cyclic_chars[i:] + self.cyclic_chars[:i]
            permutations.add(''.join(perm))
        return permutations

    def is_conjugate(self, other: 'FreeGroupWord') -> bool:
        if self.length != other.length:
            return False
        self_perms = self.get_cyclic_permutations()
        other_perms = other.get_cyclic_permutations()
        return not self_perms.isdisjoint(other_perms)

    def is_equal(self, other: 'FreeGroupWord') -> bool:
        return str(self) == str(other)

    def apply_permutation(self, perm: Dict[str, str]) -> 'FreeGroupWord':
        new_chars = []
        for c in self.reduced_chars:
            if c in perm:
                mapped = perm[c]
                if isinstance(mapped, list):
                    new_chars.extend(mapped)
                else:
                    new_chars.append(mapped)
            elif INVERSE[c] in perm:
                mapped_inv = perm[INVERSE[c]]
                if isinstance(mapped_inv, list):
                    inv_mapped = invert_word(mapped_inv)
                    new_chars.extend(inv_mapped)
                else:
                    new_chars.append(INVERSE[mapped_inv])
            else:
                new_chars.append(c)
        return FreeGroupWord(new_chars)

    def apply_dehn_twist(self, twist: Dict[str, list or str]) -> 'FreeGroupWord':
        new_chars = []
        for c in self.reduced_chars:
            if c in twist:
                mapped = twist[c]
                if isinstance(mapped, list):
                    new_chars.extend(mapped)
                else:
                    new_chars.append(mapped)
            elif INVERSE[c] in twist:
                mapped_inv = twist[INVERSE[c]]
                if isinstance(mapped_inv, list):
                    inv_mapped = invert_word(mapped_inv)
                    new_chars.extend(inv_mapped)
                else:
                    new_chars.append(INVERSE[mapped_inv])
            else:
                new_chars.append(c)
        return FreeGroupWord(new_chars)


def is_same_under_permutations(w1: FreeGroupWord, w2: FreeGroupWord) -> bool:
    for perm in PERMUTATIONS.values():
        w1_perm = w1.apply_permutation(perm)
        if w1_perm.is_equal(w2) or w1_perm.is_conjugate(w2):
            return True
    return False


def generate_first_batch_with_details(omega: FreeGroupWord) -> Tuple[List[FreeGroupWord], List[Dict]]:
    omega_len = omega.get_length()
    candidates = []
    twist_details = []
    for twist_name, twist in DEHN_TWISTS.items():
        twisted = omega.apply_dehn_twist(twist)
        twist_info = {
            'twist_name': twist_name,
            'input_word': str(omega),
            'output_raw': ''.join(twisted.raw_chars),
            'output_reduced': ''.join(twisted.reduced_chars),
            'output_word': str(twisted),
            'length': twisted.get_length(),
            'is_valid': twisted.get_length() == omega_len and not is_same_under_permutations(twisted, omega)
        }
        twist_details.append(twist_info)
        if twist_info['is_valid']:
            candidates.append(twisted)

    unique_candidates = []
    for w in candidates:
        duplicate = False
        for kept in unique_candidates:
            if is_same_under_permutations(w, kept):
                duplicate = True
                break
        if not duplicate:
            unique_candidates.append(w)
    return unique_candidates, twist_details


def generate_next_batch_with_details(prev_batch: List[FreeGroupWord], history: List[FreeGroupWord], batch_num: int) -> \
        Tuple[List[FreeGroupWord], List[List[Dict]]]:
    if not prev_batch:
        return [], []
    all_candidates = []
    batch_details = []
    for w in prev_batch:
        word_details = []
        for twist_name, twist in DEHN_TWISTS.items():
            twisted = w.apply_dehn_twist(twist)
            duplicate_with_history = False
            duplicate_with = ""
            for h in history:
                if is_same_under_permutations(twisted, h):
                    duplicate_with_history = True
                    duplicate_with = str(h)
                    break
            twist_info = {
                'twist_name': twist_name,
                'input_word': str(w),
                'output_raw': ''.join(twisted.raw_chars),
                'output_reduced': ''.join(twisted.reduced_chars),
                'output_word': str(twisted),
                'length': twisted.get_length(),
                'is_valid': twisted.get_length() == history[0].get_length() and not duplicate_with_history,
                'duplicate': duplicate_with_history,
                'duplicate_with': duplicate_with
            }
            word_details.append(twist_info)
            if twist_info['is_valid']:
                all_candidates.append(twisted)
        batch_details.append(word_details)

    unique_results = []
    for w in all_candidates:
        duplicate = False
        for kept in unique_results:
            if is_same_under_permutations(w, kept):
                duplicate = True
                break
        if not duplicate:
            unique_results.append(w)
    return unique_results, batch_details


def generate_nw1_new_words(shortest_word: str) -> List[str]:
    chars = parse_expression(shortest_word)
    omega = FreeGroupWord(chars)
    history = [omega]
    batch_results = {1: []}
    current_batch, _ = generate_first_batch_with_details(omega)
    batch_results[1] = current_batch
    if current_batch:
        history.extend(current_batch)
    batch_num = 2
    while batch_num <= 10:
        prev_batch = batch_results[batch_num - 1]
        current_batch, _ = generate_next_batch_with_details(prev_batch, history, batch_num)
        if not current_batch:
            break
        batch_results[batch_num] = current_batch
        history.extend(current_batch)
        batch_num += 1
    all_new_words = []
    for batch in batch_results.values():
        all_new_words.extend([str(w) for w in batch])
    return all_new_words


def merge_and_deduplicate(nsw_result: str, nw1_words: List[str]) -> Set[str]:
    def to_free_group_word(s: str) -> FreeGroupWord:
        chars = parse_expression(s)
        return FreeGroupWord(chars)

    merged_set = set()
    nsw_word = to_free_group_word(nsw_result)
    merged_set.add(str(nsw_word))

    for word_str in nw1_words:
        current_word = to_free_group_word(word_str)
        duplicate = False
        for existing_str in merged_set:
            existing_word = to_free_group_word(existing_str)
            if is_same_under_permutations(current_word, existing_word):
                duplicate = True
                break
        if not duplicate:
            merged_set.add(str(current_word))
    return merged_set


def apply_whitehead_permutation(word_str: str) -> Set[str]:
    perm_results = set()
    chars = parse_expression(word_str)
    base_word = FreeGroupWord(chars)

    for perm in PERMUTATIONS.values():
        perm_word = base_word.apply_permutation(perm)
        perm_results.add(str(perm_word))
    return perm_results


def get_conjugate_words(word_str: str) -> Set[str]:
    conjugate_set = set()
    chars = parse_expression(word_str)
    word = FreeGroupWord(chars)
    cyclic_perms = word.get_cyclic_permutations()

    for perm_str in cyclic_perms:
        perm_chars = list(perm_str)
        reduced_chars = reduce_word(perm_chars)
        cyclic_reduced_chars = cyclic_reduce(reduced_chars)
        conjugate_set.add(''.join(cyclic_reduced_chars))
    return conjugate_set


def get_oriented_closed_curve_words(base_set: Set[str]) -> List[str]:
    def char_key(c):
        if c == 'a':
            return 0
        elif c == 'b':
            return 1
        elif c == 'A':
            return 2
        elif c == 'B':
            return 3
        else:
            return 4

    def word_key(word):
        return [char_key(c) for c in word]

    final_set = set(base_set)

    for word in base_set:
        perm_words = apply_whitehead_permutation(word)
        final_set.update(perm_words)

        for perm_word in perm_words:
            conjugate_words = get_conjugate_words(perm_word)
            final_set.update(conjugate_words)

        original_conjugate = get_conjugate_words(word)
        final_set.update(original_conjugate)

    sorted_list = sorted(final_set, key=word_key)
    return sorted_list


# ===================== Core Functional Functions =====================
def get_oriented_closed_curve_repr(word_str: str) -> str:
    """
    Get the shortest word representation of oriented closed curve for a single word
    :param word_str: Input free group word
    :return: Shortest word representation of oriented closed curve
    """
    # 1. Get shortest cyclically reduced word
    shortest = get_shortest_cyclic_reduced_word(word_str)
    # 2. Generate nw1 new words
    nw1_words = generate_nw1_new_words(shortest)
    # 3. Merge and deduplicate
    merged = merge_and_deduplicate(shortest, nw1_words)
    # 4. Get oriented closed curve shortest word set and take first element as representation
    oriented = get_oriented_closed_curve_words(merged)
    return oriented[0] if oriented else ""


# ===================== Main Program Entry =====================
def main():
    print("===== Smallest-Shortest Word Calculation Algarithm in Automorphism Orbit of Free Group F(a,b) =====")
    print("Usage Instructions:")
    print("1. Input format: a, b, A=a^-1, B=b^-1;Input examples: abaB, abAB")
    print("2. Enter 'quit' to exit the program")
    print("========================================\n")

    while True:
        # Get user input
        user_input = input("Please enter a word in the free group: ").strip()

        # Exit logic
        if user_input.lower() == 'quit':
            print("Program exited!")
            return

        # Empty input handling
        if not user_input:
            print("Error: Input cannot be empty, please re-enter!")
            continue

        try:
            # Calculate oriented closed curve shortest word representation
            oriented_repr = get_oriented_closed_curve_repr(user_input)

            # Output result
            print(f"The smallest-shortest word in the automorphism group orbit of this word: {oriented_repr}")
            print("========================================\n")

        except ValueError as e:
            print(f"Error: {str(e)}")
            print("Please check the input format and try again!\n")
        except Exception as e:
            print(f"Program exception: {str(e)}")
            print("Please check the input format and try again!\n")


if __name__ == "__main__":
    main()