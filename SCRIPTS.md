# AstroBatch - Справочник по скриптам

Набор инструментов для пакетной обработки астрономических FITS-изображений.

---

## Краткий список скриптов

| Скрипт | Назначение |
|--------|------------|
| **add.py** | Сложение: `result = input + operand + offset` |
| **sub.py** | Вычитание: `result = input - operand + offset` |
| **mul.py** | Умножение: `result = input * operand * scale` |
| **div.py** | Деление: `result = (input / operand) * scale` |
| **arith.py** | Универсальная арифметика (add/sub/mul/div в одном скрипте) |
| **sum.py** | Суммирование стека изображений с обновлением EXPTIME |
| **med.py** | Медианное комбинирование (тайловый параллельный режим) |
| **calibrate.py** | Калибровка: dark, bias, flat, cosmetic correction |
| **autocalibrate.py** | Автокалибровка по деревьям dark/flat с подбором по EXPTIME/FILTER |
| **normalize.py** | Нормализация яркости между кадрами (линейная регрессия) |
| **ngain.py** | Нормализация умножением (приведение медианы к целевому значению) |
| **noffset.py** | Нормализация смещением (приведение медианы к целевому значению) |
| **autoflat.py** | Выравнивание градиента на флэтах (полиномиальная коррекция) |
| **cosme.py** | Коррекция горячих пикселей по списку |
| **make_cosme.py** | Генерация списка горячих пикселей из дарка |
| **makedark.py** | Создание master dark и cosme.lst из сырых дарков |
| **makeflat.py** | Создание master flat по фильтрам из сырых флэтов |
| **newflat.py** | Добавление записи в лог обслуживания (для autocalibrate --flatlog) |
| **darkopt.py** | Оптимизированное вычитание дарка (подбор коэффициента K) |
| **sortfits.py** | Сортировка FITS по времени, разбиение на сессии |
| **autosolve.py** | Астрометрическое решение (WCS), ректификация, выравнивание |
| **fft_align.py** | FFT-выравнивание кадров (поворот, масштаб, сдвиг) |

---

## Общие правила для всех скриптов

### Форматы входных данных (input_spec)

- **Одиночный файл**: `image.fit`
- **Нумерованная последовательность**: `light0001.fit` → автоматически находит light0002.fit, light0003.fit, ...
- **Маска с wildcards**: `*.fit`, `light_*.fit`
- **Файл-список**: `@list.txt` или `list.txt` (по одному пути на строку)

### Числовые константы как операнды

Арифметические скрипты (add, sub, mul, div, arith) поддерживают числовые константы в качестве операндов. Константа интерпретируется как "виртуальный файл", заполненный этим значением.

**Допустимые форматы констант** (строгий парсинг):
- Целые числа: `123`, `-456`, `+789`
- Десятичные: `123.456`, `-0.5`, `+1.0`

**НЕ допускаются**:
- Экспоненциальная нотация: `1e5`, `1E-3`
- Специальные значения: `inf`, `-inf`, `nan`
- Разделители разрядов: `1_000`, `1,000`
- Неполные числа: `.5`, `5.`

**Важно**: Хотя бы один аргумент должен быть файлом — из него берутся размеры изображения и FITS-заголовок.

**Применение**: Например может быть полезно для калибровки на термостабильных CCD, где BIAS-кадры можно заменить нулём или известной константой смещения.

### Форматы выходных данных (output_spec)

- **Одиночный файл**: `output.fit`
- **Нумерованный шаблон**: `out0001.fit` → автоинкремент для каждого входного файла

### Типы данных

**Поддерживаемые форматы:**
- Целочисленные: int8, int16, int32, int64 (signed/unsigned)
- Вещественные: float32, float64

**Выходная конвертация:**
- Целочисленные FITS: результат ограничивается диапазоном типа (clamp), округляется
- Вещественные FITS: сохраняются как есть
- NaN, Inf, -Inf: заменяются на 0 (никогда не записываются в выходные файлы)

**Когда используется float64:**
- Умножение/деление на коэффициенты ~1.0 (нормализация, flat-коррекция)
- Вычисление среднего (накопление ошибок округления)
- Операции с дробными коэффициентами

**Когда float64 НЕ нужен:**
- Сложение/вычитание целых значений или изображений
- Медиана (выбор элемента, без арифметики)
- Min/max операции

---

## Подробное описание скриптов

---

### add.py

**Назначение**: Сложение значения или изображения с входными кадрами.

**Формула**: `result = input + operand + offset`

**Синтаксис**:
```
add.py input_spec output_spec operand [offset]
```

**Параметры**:
- `input_spec` — входные файлы
- `output_spec` — выходные файлы
- `operand` — числовая константа ИЛИ FITS-файл ИЛИ нумерованный шаблон FITS
- `offset` — опциональное числовое смещение (по умолчанию 0)

**Примеры**:
```batch
add light0001.fit cal0001.fit 100
add *.fit out0001.fit bias.fit
add light0001.fit result0001.fit dark0001.fit 500
add image.fit result.fit 1024           :: добавить константу ко всем пикселям
```

---

### sub.py

**Назначение**: Вычитание значения или изображения из входных кадров.

**Формула**: `result = input - operand + offset`

**Синтаксис**:
```
sub.py input_spec output_spec operand [offset]
```

**Параметры**:
- `input_spec` — входные файлы
- `output_spec` — выходные файлы
- `operand` — числовая константа ИЛИ FITS-файл (значение для вычитания)
- `offset` — опциональное смещение (по умолчанию 0)

**Примеры**:
```batch
sub light0001.fit dark_sub0001.fit master_dark.fit
sub raw.fit calibrated.fit dark.fit 100
sub light0001.fit cal0001.fit 0         :: вычитание нуля (для термостабильных CCD без bias)
sub light0001.fit cal0001.fit 1024      :: вычитание константы смещения
```

---

### mul.py

**Назначение**: Умножение входных кадров на значение или изображение.

**Формула**: `result = input * operand * scale`

**Синтаксис**:
```
mul.py input_spec output_spec operand [scale]
```

**Параметры**:
- `input_spec` — входные файлы
- `output_spec` — выходные файлы
- `operand` — множитель (числовая константа или FITS)
- `scale` — дополнительный множитель (по умолчанию 1.0)

**Примеры**:
```batch
mul light0001.fit scaled0001.fit 2.5
mul image.fit result.fit gain_map.fit
mul light0001.fit doubled0001.fit 2     :: умножение на константу
```

---

### div.py

**Назначение**: Деление входных кадров на значение или изображение.

**Формула**: `result = (input / operand) * scale`

**Особенности**: Деление на ноль даёт 0.

**Синтаксис**:
```
div.py input_spec output_spec operand [scale]
```

**Параметры**:
- `input_spec` — входные файлы
- `output_spec` — выходные файлы
- `operand` — делитель (числовая константа или FITS)
- `scale` — множитель результата (по умолчанию 1.0)

**Примеры**:
```batch
div light0001.fit flat_div0001.fit master_flat.fit 5000
div image.fit normalized.fit 32768      :: деление на константу (нормализация к диапазону)
div light0001.fit norm0001.fit 1 5000   :: деление на 1 с масштабированием (эквивалент mul)
```

---

### arith.py

**Назначение**: Универсальный арифметический инструмент (объединяет add/sub/mul/div).

**Синтаксис**:
```
arith.py (add|sub|mul|div) input_spec output_spec operand [param]
```

**Операции**:
- `add` → `result = input + operand + offset`
- `sub` → `result = input - operand + offset`
- `mul` → `result = input * operand * scale`
- `div` → `result = (input / operand) * scale`

**Параметры**:
- `operand` — числовая константа ИЛИ FITS-файл ИЛИ нумерованный шаблон
- `param` — для add/sub это offset (по умолчанию 0), для mul/div это scale (по умолчанию 1)

**Примеры**:
```batch
arith add light0001.fit out0001.fit 100
arith sub light0001.fit cal0001.fit dark.fit 500
arith sub light0001.fit cal0001.fit 0           :: вычитание нуля (замена bias для термостабильных CCD)
arith div image.fit result.fit flat.fit 5000
arith mul light0001.fit scaled0001.fit 2.5      :: умножение на константу
```

---

### sum.py

**Назначение**: Суммирование стека FITS-изображений.

**Особенности**:
- Обновляет EXPTIME (суммарная экспозиция)
- Обновляет DATE-OBS (самое раннее время)
- Для целочисленных FITS: масштабирует результат до максимума типа
- Для вещественных FITS: вычисляет среднее

**Синтаксис**:
```
sum.py input_spec output.fit
```

**Параметры**:
- `input_spec` — входные файлы (последовательность, маска или одиночный файл)
- `output.fit` — выходной файл

**Примеры**:
```batch
sum light0001.fit stacked.fit
sum *.fit combined.fit
```

---

### med.py

**Назначение**: Медианное комбинирование стека FITS-изображений.

**Особенности**:
- Тайловый параллельный режим для больших изображений
- Точная попиксельная медиана
- Использует все ядра CPU

**Синтаксис**:
```
med.py input_spec output.fit [--tile N]
```

**Параметры**:
- `input_spec` — входные файлы
- `output.fit` — выходной файл
- `--tile N` — размер тайла в пикселях (по умолчанию 2048, 0 = без тайлов)

**Примеры**:
```batch
med dark0001.fit master_dark.fit
med flat*.fit master_flat.fit --tile 1024
med @list.txt median.fit --tile 0
```

---

### calibrate.py

**Назначение**: Калибровка изображений (вычитание dark/bias, деление на flat, cosmetic correction).

**Формула**: `result = ((raw - bias) - dark * OPTIMIZ) * K / flat`

**Синтаксис**:
```
calibrate.py input_spec output_spec -d dark.fit [-b bias.fit] [-f flat.fit K] [-c cosme.lst] [-optimize|-o]
```

**Параметры**:
- `input_spec` — сырые файлы
- `output_spec` — выходные файлы
- `-d dark.fit` — master dark (обязательно)
- `-b bias.fit` — master bias (опционально)
- `-f flat.fit K` — master flat и множитель K (опционально)
- `-c cosme.lst` — список горячих пикселей (опционально)
- `-optimize` или `-o` — оптимизация коэффициента дарка для каждого кадра

**Примеры**:
```batch
calibrate raw0001.fit cal0001.fit -d dark10s.fit
calibrate raw0001.fit cal0001.fit -d dark.fit -b bias.fit -f flat.fit 5000
calibrate raw0001.fit cal0001.fit -d dark.fit -f flat.fit 5000 -c cosme.lst -o
```

---

### autocalibrate.py

**Назначение**: Автоматическая калибровка с подбором dark/flat по EXPTIME и FILTER.

**Особенности**:
- Ищет darks по EXPTIME и JD (ближайший по времени)
- Ищет flats по FILTER и JD (ближайший, не более N дней в будущем)
- Автоматически находит cosme.lst в папке darks (по имени: dark600s.fit → cosme600s.lst)
- Множитель MULT вычисляется автоматически как медиана flat-кадра

**Синтаксис**:
```
autocalibrate.py [options] rawfiles.fit out_path dark_path flat_path
```

**Параметры**:
- `rawfiles.fit` — сырые файлы (файл, список, маска)
- `out_path` — базовое имя выходных файлов
- `dark_path` — директория с master darks и cosme*.lst
- `flat_path` — директория с master flats

**Опции**:
- `--bestflat` — автоподбор лучшего flat для каждой сессии (по минимуму σ после деления)
- `--debug` — сохранять превью в ./debug/ для визуальной проверки подбора flat
- `--flat-future-days N` — макс. дней в будущем для поиска flat (по умолчанию 2)
- `--flatlog FILE` — CSV-лог обновления флэтов для строгого подбора по интервалам

**Формула калибровки**: `result = ((raw - dark) * MULT) / flat`

Порядок: dark subtraction → cosmetic correction → flat division

**Формат выходных файлов**: `<base>_exp<EXPTIME>_<FILTER>_<N>.fit`

**Режим --bestflat**:
- Группирует файлы по фильтру и сессии (noon-to-noon по локальному времени)
- Для каждой группы создаёт превью: downscale 4×4 после dark subtraction
- Для каждого flat: медианный блур 16px → downscale → деление → удаление градиента → σ
- Выбирает flat с минимальной σ (меньше артефактов пыли/виньетирования)

**Режим --flatlog**:
- Формат CSV: `DATETIME_UTC,CAMERA_ID[,COMMENT]`
- Camera ID сопоставляется как подстрока с INSTRUME
- Логика выбора flat:
  1. Сначала ищет в текущем интервале (между записями в логе)
  2. Если в интервале несколько флэтов — берёт ближайший (не новее +N дней)
  3. Если с +N дней ничего, но в интервале есть флэт — берёт его (интервал = приоритет)
  4. Если интервал пуст — ищет глобально (по +N дней, как без flatlog)
  5. L fallback в конце, если точный фильтр не найден

**Примеры**:
```batch
:: Базовое использование
autocalibrate raw*.fit calibrated darks/ flats/

:: С автоподбором лучшего flat
autocalibrate --bestflat raw*.fit calibrated darks/ flats/

:: С отладочными превью
autocalibrate --bestflat --debug raw*.fit calibrated darks/ flats/

:: С логом обслуживания камеры
autocalibrate --flatlog maintenance.csv raw*.fit cal darks/ flats/

:: Комбинация опций
autocalibrate --bestflat --flat-future-days 3 --flatlog maint.csv raw*.fit out darks/ flats/
```

---

### newflat.py

**Назначение**: Добавление записи в лог обслуживания камеры (для --flatlog в autocalibrate).

**Применение**: Отмечает момент, когда flat-кадры становятся невалидными (чистка сенсора, замена фильтров, и т.д.)

**Синтаксис**:
```
newflat.py --camera CAMERA_ID --log FILE [--date DATETIME] [--comment TEXT]
```

**Параметры**:
- `--camera` — идентификатор камеры (должен быть подстрокой INSTRUME в FITS)
- `--log` — путь к CSV-файлу лога
- `--date` — опционально: дата/время в ISO формате (по умолчанию текущее UTC)
- `--comment` — опционально: комментарий

**Формат лога**:
```csv
# Maintenance log
DATETIME_UTC,CAMERA_ID,COMMENT
2024-05-18T14:30:00,2600MM,Cleaned sensor
2024-12-21T10:00:00,2600MM,Changed dust cover
```

**Примеры**:
```batch
:: Добавить запись с текущим временем
newflat --camera 2600MM --log maintenance.csv

:: С указанием даты
newflat --camera 2600MM --log maintenance.csv --date 2024-05-18T14:30:00

:: С комментарием
newflat --camera "ASI2600" --log maint.csv --comment "Sensor cleaning"
```

---

### normalize.py

**Назначение**: Нормализация яркости изображений относительно референсного кадра.

**Модель**: `I = B * R + C` → нормализованный результат: `(I - C) / B`

**Синтаксис**:
```
normalize.py input_spec output_spec [basefile.fit] [method]
```

**Параметры**:
- `input_spec` — входные файлы
- `output_spec` — выходные файлы
- `basefile.fit` — референсный кадр (по умолчанию первый входной)
- `method` — метод нормализации:
  - `1` — линейная регрессия (по умолчанию)
  - `2` — робастная регрессия (sigma-clipping)
  - `3` — глобальная итеративная нормализация всех кадров

**Примеры**:
```batch
normalize light0001.fit norm0001.fit
normalize light0001.fit norm0001.fit reference.fit 2
normalize *.fit normalized0001.fit 3
```

---

### ngain.py

**Назначение**: Нормализация умножением — приводит медиану каждого кадра к целевому значению.

**Формула**: `result = input * (target_median / current_median)`

**Аналог**: NGAIN в Iris.

**Особенности**:
- Поддерживает все типы данных: signed/unsigned int 8/16/32/64, float 32/64
- Все вычисления в float64
- Корректная конвертация обратно с контролем границ диапазона
- При нулевой медиане изображение копируется без изменений

**Синтаксис**:
```
ngain.py input_spec output_spec target_median
```

**Параметры**:
- `input_spec` — входные файлы
- `output_spec` — выходные файлы
- `target_median` — целевое значение медианы для всех выходных изображений

**Примеры**:
```batch
ngain flat0001.fit norm_flat0001.fit 10000
ngain *.fit normalized0001.fit 5000
ngain @list.txt out0001.fit 32768
```

---

### noffset.py

**Назначение**: Нормализация смещением — приводит медиану каждого кадра к целевому значению путём добавления константы.

**Формула**: `result = input + (target_median - current_median)`

**Аналог**: NOFFSET в Iris.

**Особенности**:
- Поддерживает все типы данных: signed/unsigned int 8/16/32/64, float 32/64
- Все вычисления в float64
- Корректная конвертация обратно с контролем границ диапазона
- Сохраняет относительные яркости пикселей (в отличие от ngain)

**Синтаксис**:
```
noffset.py input_spec output_spec target_median
```

**Параметры**:
- `input_spec` — входные файлы
- `output_spec` — выходные файлы
- `target_median` — целевое значение медианы для всех выходных изображений

**Примеры**:
```batch
noffset bias0001.fit norm_bias0001.fit 1000
noffset *.fit normalized0001.fit 5000
noffset light0001.fit out0001.fit 10000
```

---

### autoflat.py

**Назначение**: Выравнивание градиента на flat-кадрах (полиномиальная коррекция фона).

**Алгоритм**:
1. Расширение маски нулевых пикселей
2. Медианная фильтрация
3. Min-биннинг для получения фона
4. Полиномиальная аппроксимация поверхности
5. Коррекция: `result = input - model + offset`

**Синтаксис**:
```
autoflat.py [-d] input_spec output_spec [poly_order]
```

**Параметры**:
- `-d` — debug режим (сохраняет промежуточные изображения)
- `input_spec` — входные файлы
- `output_spec` — выходные файлы
- `poly_order` — порядок полинома (по умолчанию 1: плоскость)

**Примеры**:
```batch
autoflat flat0001.fit corrected0001.fit
autoflat -d flat.fit debug_flat.fit 2
```

---

### cosme.py

**Назначение**: Коррекция горячих пикселей по списку координат.

**Алгоритм**: Заменяет горячие пиксели средним значением соседей (3×3).

**Синтаксис**:
```
cosme.py input_spec output_spec cosme.lst [-t]
```

**Параметры**:
- `input_spec` — входные файлы
- `output_spec` — выходные файлы
- `cosme.lst` — список горячих пикселей (формат: `P x y`)
- `-t` — тестовый режим: создаёт маску вместо коррекции

**Формат cosme.lst**:
```
P 123 456
P 789 101
# комментарий
```

**Примеры**:
```batch
cosme light0001.fit clean0001.fit cosme.lst
cosme image.fit mask.fit cosme.lst -t
```

---

### make_cosme.py

**Назначение**: Генерация списка горячих пикселей из master dark.

**Алиас**: `find_hot` (команда `find_hot` вызывает этот же скрипт)

**Алгоритм**: Находит 10000 самых ярких пикселей в изображении.

**Синтаксис**:
```
make_cosme.py input.fit cosme.lst
```

**Параметры**:
- `input.fit` — master dark
- `cosme.lst` — выходной файл списка

**Примеры**:
```batch
make_cosme master_dark.fit cosme.lst
```

---

### makedark.py

**Назначение**: Автоматическое создание master dark и списков косметической коррекции из сырых дарков.

**Метаскрипт**: Использует sub.py, med.py, make_cosme.py.

**Алгоритм**:
1. Сканирует входные файлы, отбирает с `IMAGETYP='Dark Frame'`
2. Группирует по времени экспозиции (`EXPTIME` или `EXPOSURE`)
3. Для каждой группы:
   - Вычитает bias из каждого дарка
   - Медианно комбинирует
   - Генерирует список горячих пикселей

**Синтаксис**:
```
makedark.py input_spec [bias_spec]
```

**Параметры**:
- `input_spec` — директория ИЛИ маска ИЛИ последовательность дарков
- `bias_spec` — опционально: числовая константа ИЛИ FITS-файл ИЛИ список/маска bias-файлов (по умолчанию 0)

**Выходные файлы** (в текущую директорию):
- `dark<exp>.fit` — master dark для каждой экспозиции
- `cosme<exp>.lst` — список горячих пикселей
- `bias.fit` — master bias (если был указан список bias-файлов)

**Формат имени**: `dark300s.fit` для >= 1с, `dark500ms.fit` для < 1с

**Примеры**:
```batch
makedark C:\Darks                          :: все дарки из папки, bias=0
makedark dark0001.fit                      :: последовательность дарков, bias=0
makedark *.fit 1024                        :: маска, bias=константа
makedark C:\Darks master_bias.fit          :: с готовым master bias
makedark C:\Darks bias0001.fit             :: создаст bias.fit из последовательности
```

---

### makeflat.py

**Назначение**: Автоматическое создание master flat для каждого фильтра из сырых флэтов.

**Метаскрипт**: Использует sub.py, ngain.py, med.py, cosme.py, makedark.py.

**Алгоритм**:
1. Сканирует входные файлы, отбирает с `IMAGETYP='Flat Frame'`
2. Группирует по фильтру (`FILTER`)
3. Валидирует: все файлы в группе должны иметь одинаковую экспозицию
4. Ищет `dark<exp>.fit` и `cosme<exp>.lst` (сначала в текущей, потом во входной директории)
5. Если дарки не найдены — запускает `makedark` автоматически
6. Для каждой группы фильтра:
   - Вычитает соответствующий dark
   - Нормализует (ngain) к target_median
   - Медианно комбинирует
   - Применяет косметическую коррекцию

**Синтаксис**:
```
makeflat.py input_spec [target_median]
```

**Параметры**:
- `input_spec` — директория ИЛИ маска ИЛИ последовательность флэтов
- `target_median` — опционально: целевая медиана для нормализации (по умолчанию 5000)

**Выходные файлы** (в текущую директорию):
- `flat_<filter>.fit` — master flat для каждого фильтра

**Коды фильтров**:
| FILTER | Код файла |
|--------|-----------|
| OIII | flat_o.fit |
| SII | flat_s.fit |
| Ha | flat_h.fit |
| L | flat_l.fit |
| R | flat_r.fit |
| G | flat_g.fit |
| B | flat_b.fit |
| (другие) | flat_<полное_имя>.fit |

**Примеры**:
```batch
makeflat C:\Flats                          :: все флэты из папки, median=5000
makeflat C:\Flats 10000                    :: с другим target_median
makeflat flat*.fit                         :: маска файлов
makeflat @list.txt 8000                    :: из списка
```

---

### darkopt.py

**Назначение**: Оптимизированное вычитание дарка с подбором коэффициента K.

**Формула**: `result = input - K * dark`

**Алгоритм подбора K**:
- Разбивает изображение на тайлы 64×64
- Исключает тайлы с нулями или насыщением
- Выбирает самый тёмный тайл
- Вычисляет K методом наименьших квадратов

**Синтаксис**:
```
darkopt.py input_spec output_spec master_dark.fit [cosme.lst]
```

**Параметры**:
- `input_spec` — входные файлы
- `output_spec` — выходные файлы
- `master_dark.fit` — master dark
- `cosme.lst` — опциональный список горячих пикселей

**Примеры**:
```batch
darkopt light0001.fit opt0001.fit master_dark.fit
darkopt light0001.fit opt0001.fit dark.fit cosme.lst
```

---

### sortfits.py

**Назначение**: Сортировка FITS-файлов по времени наблюдения, разбиение на сессии.

**Синтаксис**:
```
sortfits.py input_spec output_pattern [-s|--sessions] [--gap-hours H]
```

**Параметры**:
- `input_spec` — входные файлы
- `output_pattern` — шаблон выходных имён
- `-s`, `--sessions` — режим сессий (выходные имена: `<base>_Sssss_Fffff.fit`)
- `--gap-hours H` — разрыв в часах для разделения сессий (по умолчанию 1.0)

**Примеры**:
```batch
sortfits light0001.fit sorted0001.fit
sortfits *.fit session.fit --sessions --gap-hours 2
```

---

### autosolve.py

**Назначение**: Астрометрическое решение (WCS), репроекция на касательную плоскость (TAN), субпиксельное выравнивание.

**Особенности**:
- Использует astrometry.net (solve-field) через WSL
- Репроекция в гномоническую (TAN) проекцию
- Субпиксельное выравнивание через FFT
- Рефит WCS из .corr файлов

**Синтаксис**:
```
autosolve.py [options] input_spec output_spec
```

**Основные опции**:
- `--no-solve` — пропустить астрометрическое решение
- `--rectify` — репроекция на касательную плоскость (TAN)
- `--rect-center-ra` — RA центра проекции (градусы)
- `--rect-center-dec` — Dec центра проекции (градусы)
- `--individual` — свой центр проекции для каждого файла
- `--align` — субпиксельное выравнивание
- `--ref file.fit` — референсный кадр для выравнивания
- `--scale-low`, `--scale-high` — диапазон масштаба (arcsec/pixel)
- `--radius` — радиус поиска (градусы)

**Репроекция (--rectify)**:
Преобразует изображение в гномоническую (TAN) проекцию — проекцию на плоскость, касательную к небесной сфере. По умолчанию центр берётся из первого файла; можно задать явно через `--rect-center-ra` и `--rect-center-dec`.

**Примеры**:
```batch
autosolve light0001.fit solved0001.fit
autosolve --rectify --align --ref ref.fit light0001.fit aligned0001.fit
autosolve --rectify --rect-center-ra 180.5 --rect-center-dec 45.2 *.fit out0001.fit
```

---

### fft_align.py

**Назначение**: FFT-выравнивание кадров (поворот, масштаб, субпиксельный сдвиг).

**Особенности**:
- Фазовая корреляция для поиска сдвига
- Пирамидальный поиск угла и масштаба
- Постобработка: локальная аффинная коррекция
- Параллельная обработка

**Синтаксис**:
```
fft_align.py reference input_spec output_spec [options]
```

**Основные опции**:
- `--superfine` — пирамидальный режим высокой точности
- `--post-correction` — локальная аффинная коррекция
- `--max-angle N` — максимальный угол поиска (по умолчанию 5°)
- `--scale-delta N` — диапазон масштаба ±N (по умолчанию 0.01)
- `--flux` — режим сохранения потока (линейная интерполяция)

**Примеры**:
```batch
fft_align ref.fit light0001.fit aligned0001.fit
fft_align ref.fit *.fit out0001.fit --superfine --post-correction
fft_align ref.fit light0001.fit aligned0001.fit --flux --max-angle 2
```

---

## Зависимости

**Обязательные**:
- Python 3.6+
- numpy
- astropy

**Для продвинутых функций**:
- scipy (fft_align, autoflat, autosolve)
- reproject (autosolve ректификация)
- astrometry.net (autosolve решение)

---

## Установка

1. Запустите `Commands/setup.bat` для добавления папки Commands в PATH
2. После этого скрипты доступны как команды: `add`, `sub`, `med`, `calibrate` и т.д.

```batch
Commands\setup.bat
```

---

## Тестирование

В каждой папке скрипта есть `run.bat` для тестирования:

```batch
cd Add
run.bat
```

Тестовые данные находятся в папках `Samples1/`, `Samples2/`, `1/`.
