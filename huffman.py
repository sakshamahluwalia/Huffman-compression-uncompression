"""
Code for compressing and decompressing using Huffman compression.

=====================
Assignment #2
CSC148H1 S
Raymond Truong
Saksham Ahluwalia
Muneeb Ansari
=====================
"""

from nodes import HuffmanNode, ReadNode


class Forest:
    """ A priority queue ("forest") containing HuffmanNode trees.

    === Attributes ===
    @type trees: list[(int, HuffmanNode)]
        a list of tuples, one for every tree; the first element being the
        frequency ("weight") of the tree, and the second element being the
        tree itself
    """

    def __init__(self, freq_dict):
        """ Create a new Forest self, generating the starting forest given the
        symbols and their frequencies provided in freq_dict.

        @type self: Forest
        @type freq_dict: dict(int, int)
        @rtype: None

        >>> f = Forest({31: 2, 32: 1})
        >>> len(f.trees)
        2
        """
        # each tree starts off as its weight and its symbol stored in a leaf
        self.trees = [(freq_dict[symbol], HuffmanNode(symbol, None, None))
                      for symbol in freq_dict]

    def combine_trees(self):
        """ Combine the two lowest weighted trees in self into new trees, with
        the original two trees being the two children and the sum of the
        original two weights being the new weight.

        @type self: Forest
        @rtype: None

        >>> f = Forest({31: 2, 32: 1})
        >>> f.combine_trees()
        >>> len(f.trees) == 1
        True
        >>> f.trees[0]
        (3, HuffmanNode(None, HuffmanNode(32, None, None), HuffmanNode(31, \
None, None)))
        """
        # take the two least weighted trees
        self.trees.sort()
        t1, t2 = self.trees.pop(0), self.trees.pop(0)

        t = (t1[0] + t2[0],                         # combine their weights,
             HuffmanNode(None, t1[1], t2[1]))       # combine the trees

        # add it back to the forest
        self.trees.append(t)


# ====================
# Helper functions for manipulating bytes


def get_bit(byte, bit_num):
    """ Return bit number bit_num from right in byte.

    @param int byte: a given byte
    @param int bit_num: a specific bit number within the byte
    @rtype: int

    >>> get_bit(0b00000101, 2)
    1
    >>> get_bit(0b00000101, 1)
    0
    """
    return (byte & (1 << bit_num)) >> bit_num


def byte_to_bits(byte):
    """ Return the representation of a byte as a string of bits.

    @param int byte: a given byte
    @rtype: str

    >>> byte_to_bits(14)
    '00001110'
    """
    return "".join([str(get_bit(byte, bit_num))
                    for bit_num in range(7, -1, -1)])


def bits_to_byte(bits):
    """ Return int represented by bits, padded on right.

    @param str bits: a string representation of some bits
    @rtype: int

    >>> bits_to_byte("00000101")
    5
    >>> bits_to_byte("101") == 0b10100000
    True
    """
    return sum([int(bits[pos]) << (7 - pos)
                for pos in range(len(bits))])


# ====================
# Functions for compression


def make_freq_dict(text):
    """ Return a dictionary that maps each byte in text to its frequency.

    @param bytes text: a bytes object
    @rtype: dict{int,int}

    >>> d = make_freq_dict(bytes([65, 66, 67, 66]))
    >>> d == {65: 1, 66: 2, 67: 1}
    True
    """
    # basic counter
    frequencies = {}
    for b in text:
        if b in frequencies:
            frequencies[b] += 1
        else:
            frequencies[b] = 1

    return frequencies


def huffman_tree(freq_dict):
    """ Return the root HuffmanNode of a Huffman tree corresponding
    to frequency dictionary freq_dict.

    @param dict(int,int) freq_dict: a frequency dictionary
    @rtype: HuffmanNode

    >>> freq = {2: 6, 3: 4}
    >>> t = huffman_tree(freq)
    >>> result1 = HuffmanNode(None, HuffmanNode(3), HuffmanNode(2))
    >>> result2 = HuffmanNode(None, HuffmanNode(2), HuffmanNode(3))
    >>> t == result1 or t == result2
    True
    """
    # if there's only one symbol, add another random one so that the tree
    # can be properly generated
    if len(freq_dict) == 1:
        if 255 not in freq_dict:    # just ensure that the "random" symbol
            freq_dict[255] = 0      # we're adding doesn't happen to be the
        else:                       # same one that's already in the file
            freq_dict[254] = 0

    # generate the forest
    forest = Forest(freq_dict)

    # keep combining until one is left
    while len(forest.trees) > 1:
        forest.combine_trees()

    return forest.trees[0][1]

    # === Sources ===
    # algorithm acquired from https://www.cs.duke.edu/csed/poop/huff/info/
    # supplemented by https://en.wikipedia.org/wiki/Huffman_coding
    # tree structures verified using http://huffman.ooz.ie/


def get_codes(tree):
    """ Return a dict mapping symbols from tree rooted at HuffmanNode to codes.

    @param HuffmanNode tree: a Huffman tree rooted at node 'tree'
    @rtype: dict(int,str)

    >>> tree = HuffmanNode(None, HuffmanNode(3), HuffmanNode(2))
    >>> d = get_codes(tree)
    >>> d == {3: "0", 2: "1"}
    True
    """
    codes = {}

    def get_codes_helper(node, code, path=""):
        """ Return a dict mapping symbols from tree rooted at HuffmanNode to
        codes.

        @param HuffmanNode node: a Huffman tree rooted at node 'node'
        @param str path: the path taken so far
        @param dict code: the dict of all symbols and their paths
        @rtype: dict(int,str)
        """
        # recursively follow the left child until a symbol is found
        if node.left is not None:
            get_codes_helper(node.left, code, path + "0")
        else:
            code[node.symbol] = path

        # recursively follow the right child until a symbol is found
        if node.right is not None:
            get_codes_helper(node.right, code, path + "1")
        else:
            code[node.symbol] = path

        return code

    return get_codes_helper(tree, codes)


def number_nodes(tree):
    """ Number internal nodes in tree according to postorder traversal;
    start numbering at 0.

    @param HuffmanNode tree:  a Huffman tree rooted at node 'tree'
    @rtype: NoneType

    >>> left = HuffmanNode(None, HuffmanNode(3), HuffmanNode(2))
    >>> right = HuffmanNode(None, HuffmanNode(9), HuffmanNode(10))
    >>> tree = HuffmanNode(None, left, right)
    >>> number_nodes(tree)
    >>> tree.left.number
    0
    >>> tree.right.number
    1
    >>> tree.number
    2
    """
    # number all interior nodes, starting at 0
    label_interior_nodes(tree, 0)


def avg_length(tree, freq_dict):
    """ Return the number of bits per symbol required to compress text
    made of the symbols and frequencies in freq_dict, using the Huffman tree.

    @param HuffmanNode tree: a Huffman tree rooted at node 'tree'
    @param dict(int,int) freq_dict: frequency dictionary
    @rtype: float

    >>> freq = {3: 2, 2: 7, 9: 1}
    >>> left = HuffmanNode(None, HuffmanNode(3), HuffmanNode(2))
    >>> right = HuffmanNode(9)
    >>> tree = HuffmanNode(None, left, right)
    >>> avg_length(tree, freq)
    1.9
    """
    codes = get_codes(tree)

    # total number of symbols
    symbols = sum([freq for freq in freq_dict.values()])

    # total bits in output
    bits = sum([freq_dict[s] * len(codes[s]) for s in freq_dict.keys()])

    return bits / symbols


def generate_compressed(text, codes):
    """ Return compressed form of text, using mapping in codes for each symbol.

    @param bytes text: a bytes object
    @param dict(int,str) codes: mappings from symbols to codes
    @rtype: bytes

    >>> d = {0: "0", 1: "10", 2: "11"}
    >>> text = bytes([1, 2, 1, 0])
    >>> result = generate_compressed(text, d)
    >>> [byte_to_bits(byte) for byte in result]
    ['10111000']
    >>> text = bytes([1, 2, 1, 0, 2])
    >>> result = generate_compressed(text, d)
    >>> [byte_to_bits(byte) for byte in result]
    ['10111001', '10000000']
    """
    # combine the mapping of each symbol into a string of bits
    string = ""
    for symbol in text:
        string += codes[symbol]

    # pad the string with 0s until the length is a multiple of 8
    while len(string) % 8 != 0:
        string += "0"

    # take each block of 8 bits and convert it to a byte
    compressed = [string[i:i+8] for i in range(0, len(string), 8)]

    return bytes(bits_to_byte(bits) for bits in compressed)


def tree_to_bytes(tree):
    """ Return a bytes representation of the tree rooted at tree.

    @param HuffmanNode tree: a Huffman tree rooted at node 'tree'
    @rtype: bytes

    The representation should be based on the postorder traversal of tree
    internal nodes, starting from 0.
    Precondition: tree has its nodes numbered.

    >>> tree = HuffmanNode(None, HuffmanNode(3), HuffmanNode(2))
    >>> number_nodes(tree)
    >>> list(tree_to_bytes(tree))
    [0, 3, 0, 2]
    >>> left = HuffmanNode(None, HuffmanNode(3), HuffmanNode(2))
    >>> right = HuffmanNode(5)
    >>> tree = HuffmanNode(None, left, right)
    >>> number_nodes(tree)
    >>> list(tree_to_bytes(tree))
    [0, 3, 0, 2, 1, 0, 0, 5]
    """
    # base case: None or leaf
    if tree is None or tree.is_leaf():
        return bytes([])

    # general case: interior node
    else:
        return bytes((

            # in postorder:
            list(tree_to_bytes(tree.left)) +
            list(tree_to_bytes(tree.right)) +

            # each internal node has four bytes:
            # 1) 0 if left is a leaf, 1 otherwise
            ([0] if is_leaf(tree.left) else [1]) +

            # 2) left's symbol if it's a leaf, left's number otherwise
            ([tree.left.symbol]
             if is_leaf(tree.left)
             else [tree.left.number]) +

            # 3) 0 if right is a leaf, 1 otherwise
            ([0] if is_leaf(tree.right) else [1]) +

            # 4) right's symbol if it's a leaf, right's number otherwise
            ([tree.right.symbol]
             if is_leaf(tree.right)
             else [tree.right.number])))


def num_nodes_to_bytes(tree):
    """ Return number of nodes required to represent tree (the root of a
    numbered Huffman tree).

    @param HuffmanNode tree: a Huffman tree rooted at node 'tree'
    @rtype: bytes
    """
    return bytes([tree.number + 1])


def size_to_bytes(size):
    """ Return the size as a bytes object.

    @param int size: a 32-bit integer that we want to convert to bytes
    @rtype: bytes

    >>> list(size_to_bytes(300))
    [44, 1, 0, 0]
    """
    # little-endian representation of 32-bit (4-byte)
    # int size
    return size.to_bytes(4, "little")


def compress(in_file, out_file):
    """ Compress contents of in_file and store results in out_file.

    @param str in_file: input file whose contents we want to compress
    @param str out_file: output file, where we store our compressed result
    @rtype: NoneType
    """
    with open(in_file, "rb") as f1:
        text = f1.read()
    freq = make_freq_dict(text)
    tree = huffman_tree(freq)
    codes = get_codes(tree)
    number_nodes(tree)
    print("Bits per symbol:", avg_length(tree, freq))
    result = (num_nodes_to_bytes(tree) + tree_to_bytes(tree) +
              size_to_bytes(len(text)))
    result += generate_compressed(text, codes)
    with open(out_file, "wb") as f2:
        f2.write(result)


# ====================
# Functions for decompression


def generate_tree_general(node_lst, root_index):
    """ Return the root of the Huffman tree corresponding
    to node_lst[root_index].

    The function assumes nothing about the order of the nodes in the list.

    @param list[ReadNode] node_lst: a list of ReadNode objects
    @param int root_index: index in the node list
    @rtype: HuffmanNode

    >>> lst = [ReadNode(0, 5, 0, 7), ReadNode(0, 10, 0, 12), \
    ReadNode(1, 1, 1, 0)]
    >>> generate_tree_general(lst, 2)
    HuffmanNode(None, HuffmanNode(None, HuffmanNode(10, None, None), \
HuffmanNode(12, None, None)), \
HuffmanNode(None, HuffmanNode(5, None, None), HuffmanNode(7, None, None)))
    """
    # the root of this Huffman tree is defined by the read node at root_index
    read_node, huffman_node = node_lst[root_index], HuffmanNode()

    # the left child either becomes a leaf containing the given symbol, or
    if read_node.l_type == 0:
        huffman_node.left = HuffmanNode(read_node.l_data, None, None)

    # it becomes a new Huffman tree, recursively defined by the read node at
    # the specified index
    else:
        huffman_node.left = generate_tree_general(node_lst, read_node.l_data)

    # the right child either becomes a leaf containing the given symbol, or
    if read_node.r_type == 0:
        huffman_node.right = HuffmanNode(read_node.r_data, None, None)

    # it becomes a new Huffman tree, recursively defined by the read node at
    # the specified index
    else:
        huffman_node.right = generate_tree_general(node_lst, read_node.r_data)

    return huffman_node


def generate_tree_postorder(node_lst, root_index):
    """ Return the root of the Huffman tree corresponding
    to node_lst[root_index].

    The function assumes that the list represents a tree in postorder.

    @param list[ReadNode] node_lst: a list of ReadNode objects
    @param int root_index: index in the node list
    @rtype: HuffmanNode

    >>> lst = [ReadNode(0, 5, 0, 7), ReadNode(0, 10, 0, 12), \
    ReadNode(1, 0, 1, 0)]
    >>> generate_tree_postorder(lst, 2)
    HuffmanNode(None, HuffmanNode(None, HuffmanNode(5, None, None), \
HuffmanNode(7, None, None)), \
HuffmanNode(None, HuffmanNode(10, None, None), HuffmanNode(12, None, None)))
    """
    return recreate_tree(node_lst, root_index)[0]


def generate_uncompressed(tree, text, size):
    """ Use Huffman tree to decompress size bytes from text.

    @param HuffmanNode tree: a HuffmanNode tree rooted at 'tree'
    @param bytes text: text to decompress
    @param int size: how many bytes to decompress from text.
    @rtype: bytes
    """
    out = []

    # convert the text into a path
    path = ""
    for byte in text:
        path += byte_to_bits(byte)

    read_count, i = 0, 0
    current_node = tree

    # keep reading until we have the desired number of bytes
    while read_count != size:

        # follow the path
        if path[i] == "0":
            current_node = current_node.left
        else:
            current_node = current_node.right
        i += 1

        # once the path terminates at a leaf:
        if is_leaf(current_node):
            # record the symbol stored at that leaf and start again at the root
            out.append(current_node.symbol)
            read_count += 1
            current_node = tree

    return bytes(out)


def bytes_to_nodes(buf):
    """ Return a list of ReadNodes corresponding to the bytes in buf.

    @param bytes buf: a bytes object
    @rtype: list[ReadNode]

    >>> bytes_to_nodes(bytes([0, 1, 0, 2]))
    [ReadNode(0, 1, 0, 2)]
    """
    lst = []
    for i in range(0, len(buf), 4):
        l_type = buf[i]
        l_data = buf[i+1]
        r_type = buf[i+2]
        r_data = buf[i+3]
        lst.append(ReadNode(l_type, l_data, r_type, r_data))
    return lst


def bytes_to_size(buf):
    """ Return the size corresponding to the
    given 4-byte little-endian representation.

    @param bytes buf: a bytes object
    @rtype: int

    >>> bytes_to_size(bytes([44, 1, 0, 0]))
    300
    """
    return int.from_bytes(buf, "little")


def uncompress(in_file, out_file):
    """ Uncompress contents of in_file and store results in out_file.

    @param str in_file: input file to uncompress
    @param str out_file: output file that will hold the uncompressed results
    @rtype: NoneType
    """
    with open(in_file, "rb") as f:
        num_nodes = f.read(1)[0]
        buf = f.read(num_nodes * 4)
        node_lst = bytes_to_nodes(buf)
        # use generate_tree_general or generate_tree_postorder here
        tree = generate_tree_general(node_lst, num_nodes - 1)
        size = bytes_to_size(f.read(4))
        with open(out_file, "wb") as g:
            text = f.read()
            g.write(generate_uncompressed(tree, text, size))


# ====================
# Other functions

def improve_tree(tree, freq_dict):
    """ Improve the tree as much as possible, without changing its shape,
    by swapping nodes. The improvements are with respect to freq_dict.

    @param HuffmanNode tree: Huffman tree rooted at 'tree'
    @param dict(int,int) freq_dict: frequency dictionary
    @rtype: NoneType

    >>> left = HuffmanNode(None, HuffmanNode(99), HuffmanNode(100))
    >>> right = HuffmanNode(None, HuffmanNode(101), \
    HuffmanNode(None, HuffmanNode(97), HuffmanNode(98)))
    >>> tree = HuffmanNode(None, left, right)
    >>> freq = {97: 26, 98: 23, 99: 20, 100: 16, 101: 15}
    >>> improve_tree(tree, freq)
    >>> avg_length(tree, freq)
    2.31
    """
    # generate the best case tree given the frequency dictionary, and level
    # order traverse it to find the order that the symbols *should* be in
    optimal_symbols = get_symbols(huffman_tree(freq_dict))

    # level order traverse the original tree and replace all of the symbols
    # with the properly ordered ones without changing the shape of the tree
    level_order_replace(tree, optimal_symbols)


# ====================
# Various helper functions


def is_leaf(tree):
    """ Return whether or not tree is a leaf; i.e., whether or not tree has
    zero children.

    @type tree: HuffmanNode|None
    @rtype: bool

    >>> t = HuffmanNode(3, None, None)
    >>> is_leaf(t)
    True
    >>> t = HuffmanNode(None, HuffmanNode(5, None, None), None)
    >>> is_leaf(t)
    False
    """
    return tree is not None and tree.left is None and tree.right is None


def get_symbols(tree):
    """ Return the list of all symbols in tree, by level order traversal.

    @type tree: HuffmanNode
    @rtype: list[int]

    >>> t = huffman_tree({1: 5, 2: 4, 3: 2, 4: 6, 5: 1})
    >>> get_symbols(t) == [2, 1, 4, 5, 3]
    True
    """
    symbols = []

    # queue of nodes to visit
    visit = [tree]
    while len(visit) != 0:
        next_node = visit.pop(0)

        if next_node.left is not None:
            visit.append(next_node.left)

        if next_node.right is not None:
            visit.append(next_node.right)

        # store symbols in level order
        if next_node.symbol is not None:
            symbols.append(next_node.symbol)

    return symbols


def label_interior_nodes(tree, number):
    """ Number all of the interior nodes in tree, starting at number.

    @type tree: HuffmanNode|None
    @type number: int
    @rtype: None
    """
    # base case: None or leaf
    if tree is None or is_leaf(tree):
        pass

    # general case: interior node
    else:
        # label all of the left side, starting at the given number
        label_interior_nodes(tree.left, number)

        # next, label the right side, continuing where the numbering on the
        # left ended
        number = (tree.left.number + 1
                  if tree.left is not None and tree.left.number is not None
                  else number)
        label_interior_nodes(tree.right, number)

        # next, label the current node, continuing where the numbering on the
        # right ended
        tree.number = (tree.right.number + 1
                       if tree.right is not None
                       and tree.right.number is not None
                       else number)


def recreate_tree(node_lst, root_index):
    """ Return a rebuilt HuffmanNode rooted at root_index, given a node_lst
    containing ReadNodes. Since the list is presumed to represent a tree in
    post order, also return the farthest_left_index to see how much of tree has
    been rebuilt.

    @type node_lst: list[ReadNode]
    @type root_index: int
    @rtype: (HuffmanNode, int)

    >>> lst = [ReadNode(0, 5, 0, 7), ReadNode(0, 10, 0, 12), \
    ReadNode(1, 0, 1, 0)]
    >>> recreate_tree(lst, 2)[0]
    HuffmanNode(None, HuffmanNode(None, HuffmanNode(5, None, None), \
HuffmanNode(7, None, None)), \
HuffmanNode(None, HuffmanNode(10, None, None), HuffmanNode(12, None, None)))
    """
    # the root of this Huffman tree is defined by the read node at root_index
    read_node, huffman_node = node_lst[root_index], HuffmanNode()

    # start at the root_index given and move leftwards
    farthest_left_index = root_index

    # the right child either becomes a leaf, or
    if read_node.r_type == 0:
        huffman_node.right = HuffmanNode(read_node.r_data)

    # another Huffman tree, moving leftwards in the list until the right
    # child is filled up
    else:
        right = recreate_tree(node_lst, farthest_left_index - 1)
        huffman_node.right = right[0]
        farthest_left_index = right[1]

    # the left node either becomes a leaf, or
    if read_node.l_type == 0:
        huffman_node.left = HuffmanNode(read_node.l_data)

    # another Huffman tree, further moving leftwards in the list until the left
    # child is filled up
    else:
        left = recreate_tree(node_lst, farthest_left_index - 1)
        huffman_node.left = left[0]
        farthest_left_index = left[1]

    # return the tree as well as how far left we got into the list
    return huffman_node, farthest_left_index


def level_order_replace(tree, new_symbols):
    """ Replace the symbols in tree with new_symbols, by level order traversal.

    @type tree: HuffmanNode
    @type new_symbols: list[int]
    @rtype: None

    >>> t = huffman_tree({1: 5, 2: 4, 3: 2, 4: 6, 5: 1})
    >>> get_symbols(t) == [2, 1, 4, 5, 3]
    True
    >>> level_order_replace(t, [15, 16, 17, 18, 19])
    >>> get_symbols(t) == [15, 16, 17, 18, 19]
    True
    """
    # queue of nodes to visit
    visit = [tree]
    while len(visit) != 0:
        next_node = visit.pop(0)

        if next_node.left is not None:
            visit.append(next_node.left)

        if next_node.right is not None:
            visit.append(next_node.right)

        # replace each symbol in level order with the new one that should be
        # there
        if next_node.symbol is not None:
            next_node.symbol = new_symbols.pop(0)


if __name__ == "__main__":
    import python_ta
    python_ta.check_all(config="huffman_pyta.txt")
    import doctest
    doctest.testmod()

    import time

    mode = input("Press c to compress or u to uncompress: ")
    if mode == "c":
        fname = input("File to compress: ")
        start = time.time()
        compress(fname, fname + ".huf")
        print("compressed {} in {} seconds."
              .format(fname, time.time() - start))
    elif mode == "u":
        fname = input("File to uncompress: ")
        start = time.time()
        uncompress(fname, fname + ".orig")
        print("uncompressed {} in {} seconds."
              .format(fname, time.time() - start))