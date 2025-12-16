class Solution(dict):
    # A solution is nothing but a dict of str: Any mapping

    def __str__(self):
        return "\n".join(f"{k} {v}" for k, v in self.items())
