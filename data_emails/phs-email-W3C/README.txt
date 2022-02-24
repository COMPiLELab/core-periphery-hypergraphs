Hypergraph dataset of Enron emails derived from a corpus of crawled W3C mailing lists:
https://tides.umiacs.umd.edu/webtrec/trecent/parsed_w3c_corpus.html

A hyperedge is comprised of the sender and all recipients of an email. Each line
of the file email-W3C.txt is a space-separated list of IDs of the email
addresses.

The file core-email-W3C.txt contains the IDs of the "core" nodes. These are the
email addresses with w3.org in the domain. Each hyperedge has at least one node
in the core.

The file addresses-email-W3C.txt maps the integer IDs of nodes to email
addresses.
