# Lets the CLI also run as `python -m sardine`
import sys

from .cli import main

sys.exit(main())
