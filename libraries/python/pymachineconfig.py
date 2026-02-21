import os
import json


"""
MachineConfig class which can be used to read lasercutter and 3D
printer machine and material configurations.  The config file is
cached as a JSON export of python dictionaries.
"""
class MachineConfig:

    _config = {}  # the config as read in from the config file
    _working = {} # the, possibly modified, collapsed working config

    def __init__(self, name="PythonSCAD.json"):
        try:
            self._config = self.read(name)
            #cfg = self.gen_tst_config()
            #self.write(config=cfg)
        except:
            print("Warning: config file non existant or not read.")
            print("   Generating default.")
            self._config = self.gen_tst_config()
        self._working = self.gen_working(label="default")
        return

    def gen_color_table(self):
        color_table = [
            {"label":"L00","type":"ColorTable",
             "property":{"power":1.0,"feed":1.0,"color":0x000000}},
            {"label":"L01","type":"ColorTable",
             "property":{"power":1.0,"feed":1.0,"color":0x0000FF}},
            {"label":"L02","type":"ColorTable",
             "property":{"power":1.0,"feed":1.0,"color":0xFF0000}},
            {"label":"L03","type":"ColorTable",
             "property":{"power":1.0,"feed":1.0,"color":0x00E000}},
            {"label":"L04","type":"ColorTable",
             "property":{"power":1.0,"feed":1.0,"color":0xD0D000}},
            {"label":"L05","type":"ColorTable",
             "property":{"power":1.0,"feed":1.0,"color":0xFF8000}},
            {"label":"L06","type":"ColorTable",
             "property":{"power":1.0,"feed":1.0,"color":0x00E0E0}},
            {"label":"L07","type":"ColorTable",
             "property":{"power":1.0,"feed":1.0,"color":0xFF00FF}},
            {"label":"L08","type":"ColorTable",
             "property":{"power":1.0,"feed":1.0,"color":0xB4B4B4}},
            {"label":"L09","type":"ColorTable",
             "property":{"power":1.0,"feed":1.0,"color":0x0000A0}},
            {"label":"L10","type":"ColorTable",
             "property":{"power":1.0,"feed":1.0,"color":0xA00000}},
            {"label":"L11","type":"ColorTable",
             "property":{"power":1.0,"feed":1.0,"color":0x00A000}},
            {"label":"L12","type":"ColorTable",
             "property":{"power":1.0,"feed":1.0,"color":0xA0A000}},
            {"label":"L13","type":"ColorTable",
             "property":{"power":1.0,"feed":1.0,"color":0xC08000}},
            {"label":"L14","type":"ColorTable",
             "property":{"power":1.0,"feed":1.0,"color":0x00A0FF}},
            {"label":"L15","type":"ColorTable",
             "property":{"power":1.0,"feed":1.0,"color":0xA000A0}},
            {"label":"L16","type":"ColorTable",
             "property":{"power":1.0,"feed":1.0,"color":0x808080}},
            {"label":"L17","type":"ColorTable",
             "property":{"power":1.0,"feed":1.0,"color":0x7D87B9}},
            {"label":"L18","type":"ColorTable",
             "property":{"power":1.0,"feed":1.0,"color":0xBB7784}},
            {"label":"L19","type":"ColorTable",
             "property":{"power":1.0,"feed":1.0,"color":0x4A6FE3}},
            {"label":"L20","type":"ColorTable",
             "property":{"power":1.0,"feed":1.0,"color":0xD33F6A}},
            {"label":"L21","type":"ColorTable",
             "property":{"power":1.0,"feed":1.0,"color":0x8CD78C}},
            {"label":"L22","type":"ColorTable",
             "property":{"power":1.0,"feed":1.0,"color":0xF0B98D}},
            {"label":"L23","type":"ColorTable",
             "property":{"power":1.0,"feed":1.0,"color":0xF6C4E1}},
            {"label":"L24","type":"ColorTable",
             "property":{"power":1.0,"feed":1.0,"color":0xFA9ED4}},
            {"label":"L25","type":"ColorTable",
             "property":{"power":1.0,"feed":1.0,"color":0x500A78}},
            {"label":"L26","type":"ColorTable",
             "property":{"power":1.0,"feed":1.0,"color":0xB45A00}},
            {"label":"L27","type":"ColorTable",
             "property":{"power":1.0,"feed":1.0,"color":0x004754}},
            {"label":"L28","type":"ColorTable",
             "property":{"power":1.0,"feed":1.0,"color":0x86FA88}},
            {"label":"L29","type":"ColorTable",
             "property":{"power":1.0,"feed":1.0,"color":0xFFDB66}},
            {"label":"T1","type":"ColorTable",
             "property":{"power":1.0,"feed":1.0,"color":0xF36926}},
            {"label":"T2","type":"ColorTable",
             "property":{"power":1.0,"feed":1.0,"color":0x0C96D9}},
        ]
        return color_table

    def gen_tst_config(self):
        # FIXME: this is just an expeerimental hack to get started.

        # FIXME:
        cfg = [
            # what is the default configuration?
            {"label":"default",
             "type":"default",
             "property":{"machine":"Creality-Falcon2",
                         "head":"LED-40",
                         "material":"3mm_ply_LED-40"
                         }
             },
            # different machine configurations
            {"label":"Creality-Falcon2",
             "type":"machine",
             "property":{"max_feed":25000, #(mm/min)
                         "max_width":400, #(mm)
                         "max_len":415, #(mm)
                         "has_camera":False
                         }
             },
            {"label":"XTool-S1",
             "type":"machine",
             "property":{"max_feed":36000, #(mm/min)
                         "max_width":319, #(mm)
                         "max_len":498, #(mm)
                         "has_camera":False
                         }
             },
            # different heads which can be independent of a given
            # machine
            {"label":"LED-40",
             "type":"head",
             "property":{"max_power":40.0, # (W)
                         "wavelength":455, #(nm)
                         "has_air":True,
                         "kerf": 0.075
                         },
             },
            {"label":"LED-20",
             "type":"head",
             "property":{"max_power":20.0, # (W)
                         "wavelength":455, #(nm)
                         "has_air":True,
                         "kerf": 0.075
                         },
             },
            # materials which are dependent on the head
            # characteristics. The machines assume units in mm and
            # minutes.
            {"label":"3mm_ply_LED-40",
             "type":"material",
             "property":{"thickness":3.0, # (mm)
                         "cut_power":1000, # (0.1%)
                         "cut_feed":400, # (mm/min)
                         "engrave_power":300, # (0.1%)
                         "engrave_feed":6000 # (mm/min)
                         },
             },

            {"label":"0.25in_ply_LED-40",
             "type":"material",
             "property":{"thickness":6.35,
                         "cut_power":1000,
                         "cut_feed":200,
                         "engrave_power":300,
                         "engrave_feed":6000
                         },
             },
            {"label":"0.75in_pine_LED-40",
             "type":"material",
             "property":{"thickness":19.05,
                         "cut_power":1000,
                         "cut_feed":200,
                         "engrave_power":300,
                         "engrave_feed":6000
                         },
             },
        ]

        for c in self.gen_color_table():
            cfg.append(c)

        return cfg

    def read(self, name="PythonSCAD.json"):
        name = self.configfile(name)
        with open(name, 'r', encoding='utf-8') as f:
            cfg = json.loads(f.read())
            return cfg
        
    def write(self, config=None, name="PythonSCAD.json"):
        name = self.configfile(name)

        if config is None:
            config = self._config

        jstr = json.dumps(config, indent=4)

        if jstr is not None:
            with open(name,'w') as fout:
                fout.write(jstr)
        
        return

    def set_config(self, config):
        self._config = config
        return

    def dict(self):
        return self._config

    def get_machine(self, label=None):
        default = self._config["default"]
        machine = default["machine"]
        head = default["head"]
        print("*** FIXME: default:",default)

        if label is None:
            return default
        else:
            return default[label]

    def set_working(self, config):
        self._working = config
        return

    def get_types(self):
        types = set([x["type"] for x in self._config])
        return types
    
    def get_label_by_type(self, label):
        values = set([x["label"] for x in self._config if x["type"]==label])
        return values
    
    def get_value_by_label(self, label1, label2):
        dicts = [x for x in self._config if x["label"]==label1]
        values = [x[label2] for x in dicts]
        return values

    def get_sublabel(self, label, value):
        dicts = [x for x in self._config if x["type"]==label]
        values = [x[value] for x in dicts]
        return values

    def working_config(self):
        return self._working
    
    def gen_working(self, label="default"):
        ncfg = {}

        dcfg = self.get_sublabel(label,"property")[0]

        lbls = dcfg.keys()

        for l in lbls:
            k = dcfg[l]
            tcfg = self.get_value_by_label(k,"property")
            for i in range(len(tcfg)):
                for tk in tcfg[i].keys():
                    ncfg[tk] = tcfg[i][tk]

        return ncfg

    def modify_working_config(self, label="default"):
        dcfg = self.get_value_by_label(label, "property")

        for d in dcfg:
            for l in d.keys():
                self._working[l] = d[l]
        
        return self._working

    def configfile(self, name="PythonSCAD.json"):
        name = os.path.expanduser(name)
        xdg = os.getenv("XDG_CONFIG_HOME")
        home = os.getenv("HOME")

        if '/'==name[0] or '\\'==name[0]:
            # FIXME: need to also handle 'C:' naming
            return name

        if (xdg is not None) and (os.path.exists(os.path.join(xgd,name))):
            return os.path.join(xgd,name)

        if home is not None:
            return os.path.join(home,".config","PythonSCAD",name)

        return name

    def get_property_value(self, label, tag):
        dicts = [x for x in self._config if x["label"]==label]
        values = [x["property"][tag] for x in dicts]
        return values

    def get_subproperty_value(self, label, tag1, tag2):
        dicts = [x for x in self._config if x["label"]==label]
        values = [x["property"][tag1][tag2] for x in dicts]
        return values

    def set_property_value(self, label, tag, value):
        modified = []
        for d in self._config:
            m = d
            if m["label"]==label:
                if tag in m["property"]:
                    m["property"][tag] = value
            modified.append(m)

        self._config = modified

        return

    def set_subproperty_value(self, label, tag1, tag2, value):
        modified = []
        for d in self._config:
            m = d
            if m["label"]==label:
                if tag1 in m["property"]:
                    m["property"][tag1][tag2] = value
            modified.append(m)

        self._config = modified

        return

    # The followng functions are for manipulating the color table

    def reset_colormap(self):
        """
        reset_colormap: change potentially modified labeled power,
        feed, and color associations back to their default.  The
        default color table is compatible with LightBurn's

        """
        ct = self.gen_color_table()

        for d in ct:
            self.set_property_value(d["label"], "feed",  d["property"]["feed"])
            self.set_property_value(d["label"], "power", d["property"]["power"])
            self.set_property_value(d["label"], "color", d["property"]["color"])


    def scale_value(self, label1, label2, cfg=None):
        if cfg is None:
            cfg = self._working
        val = cfg[label1] / cfg[label2]
        return val

    def color(self, tag):
        return self.get_property_value(tag, "color")[0]

    # color2str - return the working labled color as an OpenSCAD
    #   compatible string representation of the hex value starting with a
    #   '#'
    def color2str(self, tag):
        return "#{:X}".format(self.color(tag))

    # powermap - return the working labled power
    def power(self, tag):
        return self.get_property_value(tag, "power")[0]

    # feedmap - return the working labled feed
    def feed(self, tag):
        return self.get_property_value(tag, "feed")[0]

    # setpower - overwrite the working labeled power
    def set_power(self, tag, val):
        return self.set_property_value(tag, "power", val)

    # setfeed - overwrite the working labeled feed
    def set_feed(self, tag, val):
        return self.set_property_value(tag, "feed", val)

    # setcolor - overwrite the working labeled color
    def set_color(self, tag, val):
        return self.set_property_value(tag, "color", val)

    def gen_color(self, power=-1,feed=-1):
        color = 0
        power = int(power)
        feed  = int(feed)
        if power > 1000:
            print("Error: cannot represent a power factor greater than 100.0%")
        if feed > 65535:
            print("Error: cannot represent a feed rate greater than 65535mm/min")

        power = int(power/4)
        if power != -1: color |= (power << 24)
        if feed  != -1: color |= (feed << 8)

        return color

    def gen_color2str(self, power=-1,feed=-1):
        color = self.gen_color(power,feed)
        return f"#{color:08X}"

    def gen_color2hex(self, power=-1,feed=-1):
        color = self.gen_color(power,feed)
        return f"0x{color:08X}"

