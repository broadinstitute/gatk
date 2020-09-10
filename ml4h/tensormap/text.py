import re
import logging
from typing import Dict, Tuple, Callable

import numpy as np


def token_dictionary_and_text_from_file(
        text_file: str,
        remove_special_chars: bool = True,
) -> Tuple[str, Dict[str, int]]:
    texts = []
    characters = set()
    with open(text_file) as file:
        for i, line in enumerate(file.readlines()):
            cur_line = _preprocess_sentence(line, remove_special_chars)
            [characters.add(char) for char in cur_line]
            texts.append(cur_line)
            if i % 50000 == 0:
                logging.info(f'Read {i} lines from {text_file}')
    logging.info(f'Total characters: {len(characters)}')
    char2index = dict((c, i) for i, c in enumerate(sorted(list(characters))))
    index2char = dict((i, c) for i, c in enumerate(sorted(list(characters))))
    logging.info(f'char2index:\n\n {char2index}  \n\n\n\n index2char: \n\n {index2char} \n\n\n')
    return ''.join(texts), char2index


def random_text_window_tensor(
    text: str,
    window_size: int,
    one_hot: bool = True,
) -> Callable:
    def text_from_file(tm, _, dependents={}):
        tensor = np.zeros(tm.shape, dtype=np.float32)
        random_index = np.random.randint(window_size, len(text)-window_size)
        for i, c in enumerate(text[random_index:random_index+window_size]):
            if one_hot:
                tensor[i, tm.channel_map[c]] = 1.0
            else:
                tensor[i] = tm.channel_map[c]
        if tm.dependent_map is not None:
            for i, dm in enumerate(tm.dependent_map):
                start_next_window = random_index+1+i
                dependents[dm] = np.zeros(dm.shape, dtype=np.float32)
                if dm.axes() == 1 and one_hot:
                    dependents[dm][dm.channel_map[text[start_next_window]]] = 1.0
                elif dm.axes() == 2 or (not one_hot and dm.axes() == 1):
                    for j, c in enumerate(text[start_next_window:start_next_window+dm.shape[0]]):
                        if one_hot:
                            dependents[dm][j, dm.channel_map[c]] = 1.0
                        else:
                            dependents[dm][j] = dm.channel_map[c]
                else:
                    raise ValueError(f'No method to process dependent map:{dm.name} of shape {dm.shape}.')
                logging.debug(f'\nInput text: {text[random_index:random_index+window_size]}\n Dependent: {text[start_next_window:start_next_window+dm.shape[0]]}')
        return tensor
    return text_from_file


def _preprocess_sentence(sentence, remove_special_chars):
    sentence = sentence.strip()
    if remove_special_chars:
        #replacing everything with space except (a-z, A-Z, ".", "?", "!", ",")
        sentence = re.sub(r"[^a-zA-Z?.!,]+", " ", sentence)
        sentence = sentence.strip()
    return sentence
