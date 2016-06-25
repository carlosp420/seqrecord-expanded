HISTORY
=======

0.2.5 (2016-06-25)
------------------
* Upgraded requirements to Biopython 1.67
* Removed temporal overwrite of Bio.Seq class for translation of gapped sequences.

0.2.4 (2015-11-30)
------------------
* Raise custom exception if trying to translate 'N--'.

0.2.3 (2015-11-29)
------------------
* Implemented translation of gapped sequences. It assumes that gaps are defined
  as "-".

0.2.1 (2015-09-30)
------------------
* Raises error if reading frame is not specified.

0.2.0 (2015-09-29)
------------------
* Added method to class description.
* If SeqRecord does not have the `reading_frame` parameter, it will issue a warning.
  It used to issue an exception.
* Added documentation using sphinx.
* Sequences will have '?' replaced by 'N' so translations to aminoacid will work.

0.0.2 (2015-09-13)
------------------
* Added method to degenerate sequences according to Zwick et al. method: http://www.phylotools.com/ptdegenoverview.htm

0.0.1 (2015-09-13)
------------------
* First release on PyPI.
