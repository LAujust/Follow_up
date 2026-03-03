import numpy as np
import pandas as pd
from astropy.time import Time
from astropy.coordinates import SkyCoord
import astropy.units as u
import yaml


class Scheduler:
    """
    Telescope network follow-up scheduler
    """

    def __init__(self, tnow=None):
        self.tnow = tnow if tnow is not None else Time.now()

        self.telescopes = None
        self.sources = None
        self.strategy = None

        self.filtered_sources = None
        self.schedule_table = None

    # --------------------------------------------------
    # 1. Load inputs
    # --------------------------------------------------
    def load_telescopes(self, table):
        """
        table: pandas DataFrame
        """
        self.telescopes = table.copy()

        # normalize bands
        if isinstance(self.telescopes.loc[0, 'bands'], str):
            self.telescopes['bands'] = self.telescopes['bands'].apply(
                lambda x: [b.strip() for b in x.split(',')]
            )

    def load_sources(self, table):
        """
        table: pandas DataFrame
        """
        self.sources = table.copy()

        # ensure Time format
        if not isinstance(self.sources.loc[0, 'T0'], Time):
            self.sources['T0'] = self.sources['T0'].apply(Time)

    def load_strategy(self, strategy_file):
        """
        strategy_file: yaml or json
        """
        with open(strategy_file) as f:
            self.strategy = yaml.safe_load(f)

    # --------------------------------------------------
    # 2. Source filtering
    # --------------------------------------------------
    def filter_sources(self, max_age=30 * u.day):
        """
        Remove:
        - Priority == 0
        - Tnow - T0 > max_age
        """
        df = self.sources.copy()

        age = self.tnow - df['T0']
        mask = (df['Priority'] > 0) & (age < max_age)

        self.filtered_sources = df[mask].copy()
        return self.filtered_sources

    # --------------------------------------------------
    # 3. Scheduling (FRAMEWORK)
    # --------------------------------------------------
    def schedule(self):
        """
        Core scheduling logic (framework only)
        """

        if self.filtered_sources is None:
            raise RuntimeError("Run filter_sources() first")

        obs_rows = []

        # sort by priority (high first), then by age
        srcs = self.filtered_sources.copy()
        srcs['age'] = (self.tnow - srcs['T0']).to(u.day).value
        srcs = srcs.sort_values(
            by=['Priority', 'age'],
            ascending=[False, True]
        )

        for _, src in srcs.iterrows():
            P = int(src['Priority'])
            if str(P) not in self.strategy:
                continue

            strat = self.strategy[str(P)]

            # ---- placeholder scheduling logic ----
            for _, tel in self.telescopes.iterrows():
                common_bands = set(tel['bands']) & set(strat.get('bands', []))
                if not common_bands:
                    continue

                obs_rows.append({
                    'site': tel['site'],
                    'source': src['name'],
                    'ra': src['ra'],
                    'dec': src['dec'],
                    'priority': P,
                    'bands': list(common_bands),
                    'expTime': strat.get('expTime', None),
                    'expCount': strat.get('expCount', None),
                })

        self.schedule_table = pd.DataFrame(obs_rows)
        return self.schedule_table

    # --------------------------------------------------
    # 4. Output
    # --------------------------------------------------
    def generate_obs_list(self):
        """
        Return standardized observation list
        """
        if self.schedule_table is None:
            raise RuntimeError("Run schedule() first")

        return self.schedule_table

    def export(self, filename):
        """
        Save observation list
        """
        self.schedule_table.to_csv(filename, index=False)