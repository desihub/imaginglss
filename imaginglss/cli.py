from argparse import ArgumentParser as ArgumentParser, _HelpAction
from imaginglss.analysis import targetselection
from imaginglss.analysis import tycho_veto

class CLI(object):
    def __init__(self, description,
        enable_target_plugins=False,
        enable_confidence=False,
        enable_tycho_veto=False,
        ):

        ap = ArgumentParser(add_help=False, description=description)
        self.ap = ap
        if enable_target_plugins:
            self.ap.add_argument("--extra-target-definitions", action='append', default=[],
                    help="Path to additional target definitions.")

            ns, args = self.ap.parse_known_args()
            for path in ns.extra_target_definitions:
                targetselection.load(path)
        if enable_tycho_veto:
            allvetos = [i for i in dir(tycho_veto) if not str(i).startswith( '_' )]
            ap.add_argument("--use-tycho-veto", type=str, choices=allvetos, default=None, 
                help="Apply tycho veto, must run imglss-query-tycho-veto first!")
        if enable_confidence:
            ap.add_argument("--sigma-z", type=float, default=3.0, help="apply a confidence cut in z")
            ap.add_argument("--sigma-g", type=float, default=5.0, help="apply a confidence cut in g")
            ap.add_argument("--sigma-r", type=float, default=5.0, help="apply a confidence cut in r")

        ap.add_argument("--help", action=_HelpAction)
        ap.add_argument("--conf", default=None,
            help="Path to the imaginglss config file, default is from DECALS_PY_CONFIG")

    def add_target_type_argument(self, name):
        self.ap.add_argument(name, choices=[i for i in targetselection.__all__], 
                type=lambda x:getattr(targetselection, x), 
                help="Type of Target")

    def add_argument(self, *args, **kwargs):
        self.ap.add_argument(*args, **kwargs)

    def parse_args(self):
        return self.ap.parse_args()

    def prune_namespace(self, namespace):
        d = {}
        for key, value in namespace.__dict__.items():
            # FIXME: test if the key matches target_type
            if hasattr(value, 'name'):
                value = value.name
            d[key] = value
        return d
