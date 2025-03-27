from collections import deque

class TrieNode:
    def __init__(self):
        self.children = {}
        self.fail = None
        self.output = []

class AhoCorasick:
    def __init__(self):
        self.root = TrieNode()

    def add_pattern(self, pattern, index):
        node = self.root
        for char in pattern:
            if char not in node.children:
                node.children[char] = TrieNode()
            node = node.children[char]
        node.output.append(index)

    def build_automaton(self):
        queue = deque()
        for char, child in self.root.children.items():
            child.fail = self.root
            queue.append(child)
        
        while queue:
            node = queue.popleft()
            for char, child in node.children.items():
                fail_node = node.fail
                while fail_node and char not in fail_node.children:
                    fail_node = fail_node.fail
                child.fail = fail_node.children[char] if fail_node else self.root
                child.output += child.fail.output
                queue.append(child)

    def search(self, text, patterns):
        node = self.root
        results = []
        
        for i, char in enumerate(text):
            while node and char not in node.children:
                node = node.fail
            
            if node:
                node = node.children[char]
                for pattern_index in node.output:
                    results.append((pattern_index, patterns[pattern_index], i - len(patterns[pattern_index]) + 1))
            else:
                node = self.root

        return results
