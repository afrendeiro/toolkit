#!/usr/bin/env python


def test_config_has_all_required_fields():
    from ngs_toolkit import _CONFIG as local_config
    import pkgutil
    import yaml

    def _dicts_same_keys(d1, d2):
        if type(d1) != type(d2):
            return False

        for k in d1.keys():
            if k not in d2:
                return False
            else:
                if type(d1[k]) is dict:
                    return _dicts_same_keys(d1[k], d2[k])
                else:
                    return True

    file_config = (
        pkgutil.get_data("ngs_toolkit", "config/default.yaml").decode().strip()
    )
    file_config = yaml.load(file_config)

    assert _dicts_same_keys(file_config, local_config)
