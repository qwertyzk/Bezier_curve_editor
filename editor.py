# ----------- Zuzanna Kurnicka -----------
# ---Edytor kubicznych krzywych Beziera---

# Uruchamianie programu:
# python program.py

# Podstawowa instrukcja:
# - Add curve dodaje nową krzywą
# - Kliknij LPM aby dodać nowy punkt kontrolny
# - Trzymaj ctrl żeby przesuwać punkty niezależnie od siebie
# - Trzymaj alt żeby wygładzić złączenie i przesuwać punkt wraz z jego sąsiadami

import tkinter as tk
from tkinter import filedialog
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import json


class InteractiveParametricBezierPlot:
    def __init__(self):
        self.fig, self.ax = plt.subplots(figsize=(10, 10))

        # Dodaj grid kwadratowy
        self.ax.set_aspect("equal")
        self.ax.grid(True, linestyle="--", linewidth=0.5)

        # Ustaw zakres na osiach
        self.ax.set_xlim(0, 10)
        self.ax.set_ylim(0, 10)

        # Ustaw kroki na osiach
        self.ax.set_xticks(np.arange(0, 10.1, 0.5))
        self.ax.set_yticks(np.arange(0, 10.1, 0.5))

        # Zmienne do przechowywania wielu zestawów punktów i krzywych
        self.curves = []  # Lista krzywych jako zestawów współrzędnych
        self.current_curve = -1

        # Maksymalna odległość kursora od punktu
        self.threshold = 0.1

        # Przyciski
        button_width, button_height = 0.1, 0.05
        start_pos_x, gap = 0.05, 0.03
        buttons_number = 7

        button_positions = [
            [
                start_pos_x + i * gap + i * button_width,
                0.01,
                button_width,
                button_height,
            ]
            for i in range(buttons_number)
        ]

        self.add_curve_button = plt.Button(plt.axes(button_positions[0]), "Add\ncurve")
        self.add_curve_button.on_clicked(self.add_curve)

        self.next_curve_button = plt.Button(plt.axes(button_positions[1]), "Next\ncurve")
        self.next_curve_button.on_clicked(self.next_curve)

        self.save_button = plt.Button(plt.axes(button_positions[2]), "Save\ndata")
        self.save_button.on_clicked(self.save_to_file)

        self.open_button = plt.Button(plt.axes(button_positions[3]), "Open\ndata")
        self.open_button.on_clicked(self.load_from_file)

        self.set_background_button = plt.Button(plt.axes(button_positions[4]), "Set\nbackground")
        self.set_background_button.on_clicked(self.set_background)

        self.toggle_view_button = plt.Button(plt.axes(button_positions[5]), "Toggle\nview")
        self.toggle_view_button.on_clicked(self.toggle_view)

        self.save_image_button = plt.Button(plt.axes(button_positions[6]), "Save\nimage")
        self.save_image_button.on_clicked(self.save_image)

        # Flaga określająca, czy punkty i odcinki są widoczne
        self.show_points_and_lines = True

        # Punkt aktualnie przesuwany (jego indeks w tablicy)
        self.dragging_point = None

        # Akcje
        self.cid_press = self.fig.canvas.mpl_connect(
            "button_press_event", self.on_click
        )
        self.cid_motion = self.fig.canvas.mpl_connect(
            "motion_notify_event", self.on_motion
        )
        self.cid_release = self.fig.canvas.mpl_connect(
            "button_release_event", self.on_release
        )

    # Zapisuje krzywą w postaci obrazu png
    def save_image(self, event):
        root = tk.Tk()
        root.withdraw()
        filename = filedialog.asksaveasfilename(
            defaultextension=".png", filetypes=[("PNG files", "*.png")]
        )
        if filename:
            self.show_points_and_lines = False  # Wyłącz widoczność punktów kontrolnych i pomocniczych odcinków
            self.update_plot()
            self.ax.set_axis_off()  # Wyłącz osie przed zapisaniem obrazu
            
            # Zapamiętaj tło przed jego usunięciem
            current_background = self.ax.get_images()
            # Usuń tło
            for bg in current_background:
                bg.remove()
            
            # Zapisuje obraz obszaru rysowania niezaleznie od rozmiaru okna
            extent = self.ax.get_window_extent().transformed(
                self.fig.dpi_scale_trans.inverted()
            )
            self.fig.savefig(
                filename,
                format="png",
                pad_inches=0,
                bbox_inches = extent,
                transparent=True,
            )

            # Przywróć tło
            for bg in current_background:
                self.ax.add_image(bg)
            self.show_points_and_lines = True  # Przywróć widoczność punktów kontrolnych i pomocniczych odcinków
            self.ax.set_axis_on()  # Przywróć osie
            self.update_plot()

    def save_to_file(self, event):
        root = tk.Tk()
        root.withdraw()
        filename = filedialog.asksaveasfilename(
            defaultextension=".json", filetypes=[("JSON files", "*.json")]
        )
        if filename:
            with open(filename, "w") as file:
                json.dump(self.curves, file)

    def load_from_file(self, event):
        root = tk.Tk()
        root.withdraw()
        filename = filedialog.askopenfilename(filetypes=[("JSON files", "*.json")])
        if filename:
            with open(filename, "r") as file:
                self.curves = json.load(file)
                self.current_curve = len(self.curves) - 1 if self.curves else -1
                self.update_plot()

    def set_background(self, event):
        root = tk.Tk()
        root.withdraw()
        filename = filedialog.askopenfilename(
            filetypes=[("Image files", "*.jpg *.png *.jpeg")]
        )
        if filename:
            img = mpimg.imread(filename)

            # Obliczanie proporcji obrazka
            height, width = img.shape[:2]
            aspect_ratio = width / height

            # Skalowanie i centrowanie obrazka
            if aspect_ratio > 1:  # Szerokość jest dominująca
                scaled_height = 10 / aspect_ratio
                vertical_margin = (10 - scaled_height) / 2
                extent = [0, 10, vertical_margin, 10 - vertical_margin]
            else:  # Wysokość jest dominująca
                scaled_width = 10 * aspect_ratio
                horizontal_margin = (10 - scaled_width) / 2
                extent = [horizontal_margin, 10 - horizontal_margin, 0, 10]

            self.ax.imshow(img, aspect="equal", extent=extent, alpha=0.5)
            self.fig.canvas.draw()

    def toggle_view(self, event):
        # Obsługa przycisku przełączania widoku
        self.show_points_and_lines = not self.show_points_and_lines
        self.update_plot()

    # Dodaj krzywą
    def add_curve(self, event):
        new_curve_data = {"x": [], "y": []}

        self.curves.append(new_curve_data)
        self.current_curve = len(self.curves) - 1
        self.update_plot()

    # Przełącz na następną krzywą
    def next_curve(self, event):
        if self.curves:
            self.current_curve = (self.current_curve + 1) % len(self.curves)
            self.update_plot()

    # Sprawdź, czy krzywa jest zamknięta
    def is_curve_closed(self, curve_data):
        if len(curve_data["x"]) >= 2:
            return (
                curve_data["x"][0] == curve_data["x"][-1]
                and curve_data["y"][0] == curve_data["y"][-1]
            )
        else:
            return False

    def delete_point(self, x, y, data):
        index = self.find_closest_point(x, y, data)

        if (index is not None and index % 3 == 0):  # Upewniamy się, że punkt jest punktem kontrolnym
            if index == 0 and len(data["x"]) > 2:  # Pierwszy punkt kontrolny
                del_range = slice(
                    index, index + 3
                )  # Usuwamy punkt kontrolny i dwa punkty pomocnicze
            elif (
                index == len(data["x"]) - 1 and len(data["x"]) > 2
            ):  # Ostatni punkt kontrolny
                del_range = slice(
                    index - 2, index + 1
                )  # Usuwamy punkt kontrolny i dwa punkty pomocnicze
            else:  # Wewnętrzne punkty kontrolne
                del_range = slice(max(index - 1, 0), min(index + 2, len(data["x"])))

            del data["x"][del_range]
            del data["y"][del_range]
            self.update_plot()  # Aktualizacja wykresu po usunięciu punktu

    def on_click(self, event):
        if event.inaxes != self.ax or self.current_curve == -1:
            return

        x, y = event.xdata, event.ydata
        current_data = self.curves[self.current_curve]

        if event.button == 1 and (event.key == "control" or event.key == "alt"):
            self.dragging_point = self.find_closest_point(x, y, current_data)
        elif event.button == 1 and event.key == "shift":
            self.delete_point(x, y, current_data)
        elif event.button == 1:
            if not self.is_curve_closed(
                current_data
            ):  # Sprawdź, czy krzywa jest zamknięta
                current_data["x"].append(x)
                current_data["y"].append(y)

                closest_point = self.find_closest_point(x, y, current_data)

                if closest_point is not None and closest_point == 0:
                    # Jeśli najbliższy punkt jest równy pierwszemu punktowi, zamknij krzywą
                    self.close_curve(current_data)

                if len(current_data["x"]) > 1:
                    # Dodaj dwa dodatkowe punkty po lewej i po prawej od nowo dodanego głównego punktu
                    idx = len(current_data["x"]) - 2
                    x_middle = (current_data["x"][idx] + current_data["x"][idx + 1]) / 2
                    y_middle = (current_data["y"][idx] + current_data["y"][idx + 1]) / 2

                    current_data["x"].insert(
                        -1, (2 * current_data["x"][idx] + x_middle) / 3
                    )
                    current_data["y"].insert(
                        -1, (2 * current_data["y"][idx] + y_middle) / 3
                    )

                    current_data["x"].insert(
                        -1, (2 * current_data["x"][idx + 1] + x_middle) / 3
                    )
                    current_data["y"].insert(
                        -1, (2 * current_data["y"][idx + 1] + y_middle) / 3
                    )

        self.update_plot()

    def on_release(self, event):
        if event.button == 1:
            self.dragging_point = None

    def on_motion(self, event):
        if self.dragging_point is not None and self.current_curve != -1:
            current_data = self.curves[self.current_curve]
            x, y = event.xdata, event.ydata

            if event.xdata is None or event.ydata is None:
                # Przekształcenie współrzędnych pikseli na współrzędne danych
                ax = self.fig.axes[0]
                x = ax.transData.inverted().transform([event.x, event.y])[0]
                y = ax.transData.inverted().transform([event.x, event.y])[1]
                x = np.clip(x, 0, 10)
                y = np.clip(y, 0, 10)

            if event.key == "alt":
                if np.sqrt((current_data['x'][0] - current_data['x'][-1])**2 + (current_data['y'][0] - current_data['y'][-1])**2) < self.threshold:
                    self.close_curve(current_data)


                # Logika dla wygładzania krzywej

                if self.dragging_point % 3 == 1:  # Sąsiad po lewej stronie
                    connection_point = self.dragging_point - 1

                    if connection_point - 1 >= 0:
                        right_neighbor = connection_point - 1
                    elif self.is_curve_closed(current_data):
                        right_neighbor = -2
                    else:
                        right_neighbor = None

                    if right_neighbor:
                        current_data["x"][right_neighbor] = (
                            2 * current_data["x"][connection_point] - x
                        )
                        current_data["y"][right_neighbor] = (
                            2 * current_data["y"][connection_point] - y
                        )

                elif self.dragging_point % 3 == 2:  # Sąsiad po prawej stronie
                    connection_point = self.dragging_point + 1

                    if connection_point + 1 < len(current_data["x"]):
                        left_neighbor = connection_point + 1
                    elif self.is_curve_closed(current_data):
                        left_neighbor = 1
                    else:
                        left_neighbor = None

                    if left_neighbor:
                        current_data["x"][left_neighbor] = (
                            2 * current_data["x"][connection_point] - x
                        )
                        current_data["y"][left_neighbor] = (
                            2 * current_data["y"][connection_point] - y
                        )

                elif self.dragging_point % 3 == 0:
                    connection_point = self.dragging_point
                    left_neighbor = (
                        connection_point + 1
                        if connection_point + 1 < len(current_data["x"])
                        else None
                    )
                    right_neighbor = (
                        connection_point - 1 if connection_point - 1 >= 0 else None
                    )

                    vector_x = x - current_data["x"][connection_point]
                    vector_y = y - current_data["y"][connection_point]

                    if left_neighbor is not None:
                        current_data["x"][left_neighbor] += vector_x
                        current_data["y"][left_neighbor] += vector_y

                    if right_neighbor is not None:
                        current_data["x"][right_neighbor] += vector_x
                        current_data["y"][right_neighbor] += vector_y

                    if self.is_curve_closed(current_data) and (
                        left_neighbor is None or right_neighbor is None
                    ):
                        if connection_point == len(current_data["x"]) - 1:
                            current_data["x"][0] = x
                            current_data["y"][0] = y

                            current_data["x"][1] += vector_x
                            current_data["y"][1] += vector_y
                        else:
                            current_data["x"][-1] = x
                            current_data["y"][-1] = y

                            current_data["x"][-2] += vector_x
                            current_data["y"][-2] += vector_y

            # Aktualizacja przesuwanego punktu
            current_data["x"][self.dragging_point] = x
            current_data["y"][self.dragging_point] = y

            self.update_plot()

    def close_curve(self, data):
        # Zamknij krzywą, ustawiając ostatni punkt na współrzędne pierwszego punktu
        if len(data["x"]) > 3:
            data["x"][-1] = data["x"][0]
            data["y"][-1] = data["y"][0]
            self.update_plot()

    def find_closest_point(self, x, y, data):
        # Znajdź punkt, który mieści się w ustalonej odgórnie odległości
        if not data["x"]:
            return None

        distances = np.sqrt(
            (np.array(data["x"]) - x) ** 2 + (np.array(data["y"]) - y) ** 2
        )
        close_points = np.where(distances < self.threshold)[0]

        if len(close_points) > 0:
            return close_points[0]
        else:
            return None



    def update_plot(self):
        # Usuń stare krzywe, punkty kontrolne i odcinki z osi
        for line in self.ax.get_lines():
            line.remove()

        for collection in self.ax.collections:
            collection.remove()

        # Rysuj tylko krzywe
        for curve_index, curve_data in enumerate(self.curves):
            n = len(curve_data["x"])

            if n >= 4:
                for i in range(0, n - 3, 3):
                    segment_points = list(
                        zip(curve_data["x"][i : i + 4], curve_data["y"][i : i + 4])
                    )
                    t_values = np.linspace(0, 1, 100)
                    points_interp = [
                        self.bezier_curve_function(t, segment_points) for t in t_values
                    ]
                    x_interp, y_interp = zip(*points_interp)

                    line_color = (
                        "black"
                        if not self.show_points_and_lines
                        else ("green" if curve_index == self.current_curve else "gray")
                    )
                    self.ax.plot(x_interp, y_interp, color=line_color, linewidth=2)

        if self.show_points_and_lines:
            # Rysuj punkty kontrolne i odcinki łączące
            for curve_index, curve_data in enumerate(self.curves):
                n = len(curve_data["x"])
                point_sizes = [
                    15 if i % 3 == 0 else 10 for i in range(len(curve_data["x"]))
                ]
                points_color = [
                    "gray"
                    if curve_index != self.current_curve
                    else ("red" if i % 3 == 0 else "blue")
                    for i in range(len(curve_data["x"]))
                ]
                self.ax.scatter(
                    curve_data["x"], curve_data["y"], c=points_color, s=point_sizes
                )

                if curve_index == self.current_curve:
                    color = "blue"
                else:
                    color = "gray"

                if n > 3:
                    for j in range(0, n, 3):
                        if j == 0:
                            self.ax.plot(
                                [curve_data["x"][0], curve_data["x"][1]],
                                [curve_data["y"][0], curve_data["y"][1]],
                                color=color,
                                linewidth=1,
                            )

                        elif j == n - 1:
                            self.ax.plot(
                                [curve_data["x"][n - 2], curve_data["x"][n - 1]],
                                [curve_data["y"][n - 2], curve_data["y"][n - 1]],
                                color=color,
                                linewidth=1,
                            )

                        else:
                            self.ax.plot(
                                [
                                    curve_data["x"][j - 1],
                                    curve_data["x"][j],
                                    curve_data["x"][j + 1],
                                ],
                                [
                                    curve_data["y"][j - 1],
                                    curve_data["y"][j],
                                    curve_data["y"][j + 1],
                                ],
                                color=color,
                                linewidth=1,
                            )

        self.fig.canvas.draw()

    def bezier_curve_function(self, t, control_points):
        px = [p[0] for p in control_points]
        py = [p[1] for p in control_points]

        x = (
            (1 - t) ** 3 * px[0]
            + 3 * (1 - t) ** 2 * t * px[1]
            + 3 * (1 - t) * t**2 * px[2]
            + t**3 * px[3]
        )
        y = (
            (1 - t) ** 3 * py[0]
            + 3 * (1 - t) ** 2 * t * py[1]
            + 3 * (1 - t) * t**2 * py[2]
            + t**3 * py[3]
        )

        return x, y

    def show_plot(self):
        plt.show()


interactive_parametric_bezier_plot = InteractiveParametricBezierPlot()
interactive_parametric_bezier_plot.show_plot()
