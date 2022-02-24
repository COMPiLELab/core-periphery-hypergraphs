Hypergraph dataset of Enron emails derived from William Cohen's dataset at:
https://www.cs.cmu.edu/~enron/

A hyperedge is comprised of the sender and all recipients of an email. Each line
of the file email-Enron.txt is a space-separated list of IDs of the email
addresses.

The file core-email-Enron.txt contains the IDs of the "core" nodes. These
correspond to the individuals whose email inboxes were released as part of the
investigation by the Federal Energy Regulatory Commission. Each hyperedge has at
least one node in the core.

The file addresses-email-Enron.txt maps the integer IDs of core nodes to email
addresses.
