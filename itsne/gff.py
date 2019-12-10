import gzip
import re


class GFF:

  def __init__(self,
               filename,
               feature_type=None,
               filters=None,
               gene_blacklist=None):
    self.gene_blacklist = None
    if gene_blacklist:
      self.gene_blacklist = set(
          [item.strip() for item in open(gene_blacklist, 'r').readlines()])
    self._handle = gzip.open(filename, "r")
    self.feature_type = feature_type
    self.filters = filters
    self.entries = []
    self.attr_regexes = [r"(\S+)=(\S+)", r"(\S+) \"(\S+)\""]

    while True:
      line = self._handle.readline().decode("utf-8")
      if not line:
        break

      if line.startswith("#"):
        continue

      [seqname, source, feature, start, end, score, strand, frame,
       attribute] = line.split("\t")

      if self.feature_type and feature != self.feature_type:
        continue

      entry_passes_gene_blacklist = True
      if self.gene_blacklist:
        for bad_gene in self.gene_blacklist:
          if bad_gene in attribute:
            entry_passes_gene_blacklist = False
            break

      if not entry_passes_gene_blacklist:
        continue

      result = {
          "seqname": seqname,
          "source": source,
          "feature": feature,
          "start": int(start),
          "end": int(end),
          "score": score,
          "strand": strand,
          "frame": frame
      }

      for attr_raw in [s.strip() for s in attribute.split(";")]:
        if not attr_raw or attr_raw == "":
          continue

        for regex in self.attr_regexes:
          match = re.match(regex, attr_raw)
          if match:
            key, value = match.group(1), match.group(2)
            result["attr_" + key] = value.strip()

      entry_passes_filters = True
      if entry_passes_filters and self.filters:
        for (k, v) in self.filters.items():
          if k not in result or result[k] != v:
            entry_passes_filters = False
            break

      if entry_passes_filters:
        self.entries.append(result)

  def filter(self, func):
    self.entries = filter(func, self.entries)

  def __iter__(self):
    self.i = 0
    return self

  def __next__(self):
    if self.i < len(self.entries):
      self.i += 1
      return self.entries[self.i - 1]
    else:
      raise StopIteration
