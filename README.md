# ZK_STARKS Research


![Stars](https://img.shields.io/github/stars/encrypted-systems/zk_starks?style=social)
![License](https://img.shields.io/github/license/encrypted-systems/zk_starks)
![Python](https://img.shields.io/badge/language-python-brightgreen)
![Blockchain](https://img.shields.io/badge/tech-blockchain-blue)
![Cryptography](https://img.shields.io/badge/tech-cryptography-blue)
![Requirements](https://img.shields.io/badge/dependencies-up%20to%20date-brightgreen)
![ZKP](https://img.shields.io/badge/zk-proofs-important)
![STARKs](https://img.shields.io/badge/STARKs-scaling-blue)
![Quantum Resistance](https://img.shields.io/badge/security-quantum%20resistant-brightgreen)
![Blockchain](https://img.shields.io/badge/blockchain-scalable-blue)
![Decentralized](https://img.shields.io/badge/network-decentralized-lightgreen)
![Support](https://img.shields.io/badge/support-give%20a%20star-brightgreen)
![Contributor](https://img.shields.io/badge/contributors-welcome-orange)



Автор: **Муравьев Герман Валерьевич**

***ВСЕ ПРАВА ЗАЩИЩЕНЫ ПО СРЕДСТВУ ИСПОЛЬЗОВАНИЯ MIT ЛИЦЕНЗИИ***

**Примечание:** Статья подготовлена исключительно в образовательных и научных целях.
Автор не несет ответственности за использование представленного алгоритма в целях,
связанных с незаконной деятельностью, включая, но не ограничиваясь, передачей средств,
имеющих сомнительное происхождение, через интернет-сети.
Все риски, связанные с нарушением законодательства, возлагаются на пользователя.

## ⭐ Поддержите проект
Если вам понравился этот репозиторий, пожалуйста, поставьте ⭐! Это мотивирует работать над новыми улучшениями и поддерживать проект.

---

## Оглавление
- [Введение](#введение)
- [ZK-STARKs: Основы и Руководство](#zk-starks-основы-и-руководство)
- [Глава 1: Введение в концепцию](#глава-1-введение-в-концепцию)
  - [Что такое нулевое знание?](#что-такое-нулевое-знание)
  - [Основная идея](#основная-идея)
  - [Основные свойства](#основные-свойства)
- [Глава 2: Разбор нулевого знания](#глава-2-разбор-нулевого-знания)
  - [Основная гипотеза и вычислительная задача](#основная-гипотеза-и-вычислительная-задача)
- [Глава 3: Основной принцип работы](#глава-3-основной-принцип-работы)
  - [1. Трасса вычислений и расширение низкой степени](#1-трасса-вычислений-и-расширение-низкой-степени)
  - [2. Ограничения](#2-ограничения)
  - [3. FRI Ограничения](#3-fri-ограничения)
  - [4. Фаза запроса](#4-фаза-запроса)
- [Заключение](#заключение)

---

## Введение

Меня зовут Герман, и я увлечен исследованием и созданием будущего в сфере цифровых финансов. Моя работа и интересы лежат на пересечении технологий, экономики, инвестиций и геополитики.

Я являюсь младшим криптографом и разработчиком, специализирующимся на цифровых финансах, блокчейне и вопросах безопасности в этой динамично развивающейся области.

В своих проектах я сосредоточен на инновациях и защите данных, стремясь объединить технологии и стратегическое видение для создания устойчивых решений в мире цифровой экономики.

В этой статье я намерен исследовать истоки анонимных и теневых транзакций в блокчейне, раскрывая их происхождение и влияние на экосистему.

---

## ZK-STARKs: Основы и Руководство

ZK-STARKs (Zero-Knowledge Scalable Transparent ARguments of Knowledge) — это передовая криптографическая технология, позволяющая доказать корректность вычислений без раскрытия самих данных. Этот документ поможет вам понять основные концепции ZK-STARKs, их свойства и механизм работы.

---

## Глава 1: Введение в концепцию

### Что такое нулевое знание?

Концепция нулевого знания заключается в том, чтобы одна сторона (доказатель) могла убедить другую сторону (верификатора) в истинности утверждения, **не раскрывая никаких дополнительных данных, кроме факта истинности**.

Пример: доказатель утверждает, что знает решение уравнения, но не раскрывает само решение. Верификатор убеждается в этом без получения дополнительных сведений.

Ключевые аспекты нулевого знания:
- **Проверяемость**: Верификатор убеждается в правильности утверждения.
- **Приватность**: Никакая информация, кроме истинности утверждения, не раскрывается.
- **Завершённость**: Если доказатель действует честно, то доказательство всегда принимается.

### Основная идея

ZK-STARKs предлагают масштабируемую и прозрачную систему доказательств, где:
1. **Масштабируемость**: Размер доказательства и время проверки значительно меньше, чем объем данных или сложность вычисления.
2. **Прозрачность**: Исключение необходимости доверенных сторон для генерации доказательств. Все основывается на общедоступных случайных данных (public randomness), таких как хеш-функции.

ZK-STARKs используются для доказательства корректности сложных вычислений, например:
- Проверка транзакций в блокчейнах.
- Проверка вычислений в децентрализованных системах.
- Хранение данных в приватных вычислениях.

### Основные свойства

ZK-STARKs выделяются следующими свойствами:

1. **Нулевое знание**: Доказатель убеждает верификатора в истинности вычислений, не раскрывая входных данных.
2. **Масштабируемость**: Поддерживает доказательства для больших вычислений за логарифмическое время проверки.
3. **Прозрачность**: Исключает необходимость доверенных третьих лиц благодаря использованию криптографически безопасных случайных чисел.
4. **Квантовая стойкость**: ZK-STARKs защищены от атак квантовых компьютеров, поскольку основываются на хеш-функциях вместо более уязвимых алгоритмов, таких как RSA или эллиптические кривые.
5. **Децентрализация**: ZK-STARKs идеально подходят для использования в блокчейнах и других децентрализованных системах.

---

## Глава 2: Разбор нулевого знания

В этой главе мы подробно разберем базовые идеи, лежащие в основе систем нулевого знания, включая гипотезы, которые обеспечивают их безопасность, и формулировку вычислительной задачи.

### Основная гипотеза и вычислительная задача

#### Основная гипотеза

Системы нулевого знания основаны на двух фундаментальных криптографических гипотезах:

1. **Сложность вычислительных задач**:
   Существует класс задач, которые сложно решить за приемлемое время, но их решение легко проверить.
   Пример: факторизация больших чисел, дискретный логарифм или нахождение корней полиномиальных уравнений в конечных полях.

2. **Свойства криптографических хеш-функций**:
   Для построения доказательств ZK-STARKs используются хеш-функции, которые обладают следующими свойствами:
   - **Однонаправленность**: Невозможно эффективно восстановить исходные данные из хеш-значения.
   - **Устойчивость к коллизиям**: Невозможно эффективно найти два различных входа с одинаковым хеш-значением.
   - **Псевдослучайность**: Выходные данные функции кажутся случайными и равномерно распределены.

Эти гипотезы гарантируют, что доказательства невозможно подделать или обойти.

#### Вычислительная задача

**Цель ZK-STARKs** — доказать корректность вычисления некоторой функции $f$ над данными $x$ без раскрытия $x$.

Формально:
Доказатель хочет убедить верификатора в том, что он знает такие данные $x$, что:

$$
y = f(x),
$$

где:
- $y$ — известный результат,
- $f$ — детерминированная функция.

##### Пример задачи:
- Доказатель утверждает, что он знает секретный ввод $x$, который при хешировании дает значение $y$:
  $
  y = \text{hash}(x).
  $
  Однако доказатель не раскрывает $x$.

---

## Глава 3: Основной принцип работы

ZK-STARKs основываются на последовательности этапов, которые разбиты на четыре ключевые части. Каждая из них играет важную роль в обеспечении безопасности, прозрачности и эффективности системы доказательств.


### 1. Трасса вычислений и расширение низкой степени

**[обзор данной части (ссылка)](part_1.ipynb)**

**Трасса вычислений**:
Входное вычисление разбивается на последовательность шагов, называемых *трассой вычислений*. Каждый шаг записывается в виде набора состояний, которые отражают текущий этап вычисления функции $f(x)$. Эти состояния кодируются алгебраическими выражениями.

**Расширение низкой степени**:
Чтобы защититься от манипуляций и сжать данные для проверки, состояния трассы преобразуются в полиномы низкой степени. Это достигается с помощью **Reed-Solomon кодов**:
- Трасса представляется как набор значений полинома $P(x)$ над конечным полем.
- Полином $P(x)$ расширяется на больший домен (с помощью увеличения степени), чтобы проверить, что его поведение соответствует ограничениям трассы.

Таким образом, сложное вычисление сводится к проверке свойств полинома.

### 2. Ограничения
**[обзор данной части (ссылка)](part_2.ipynb)**

На этом этапе задаются **ограничения**, которые должны удовлетворять значения трассы. Эти ограничения выражаются через полиномиальные равенства.

Пример: Если функция $f(x)$ реализована через итеративный алгоритм, каждое следующее состояние трассы можно выразить через предыдущее с использованием полиномов.

Эти полиномиальные ограничения проверяются на расширенном домене, чтобы убедиться, что трасса корректна.


### 3. FRI Ограничения

**[обзор данной части (ссылка)](part_3.ipynb)**

FRI (Fast Reed-Solomon Interactive Oracle Proof of Proximity) используется для проверки, что полином действительно имеет низкую степень без раскрытия его явного вида.

**Принцип работы FRI**:
1. Полином $P(x)$ проверяется на соответствие низкой степени.
2. Используются Merkle-деревья для обеспечения прозрачности и надежной проверки.
3. Процесс проходит через несколько раундов (слоев FRI), где каждый следующий слой содержит меньший объем данных, упрощая проверку.
На каждом этапе проверяется, что значения полинома соответствуют указанным ограничениям, постепенно уменьшая степень полинома.

Этот процесс делает систему масштабируемой и вычислительно эффективной.


#### 4. Фаза запроса

**[обзор данной части (ссылка)](part_4.ipynb)**

На завершающем этапе верификатор выполняет случайные запросы к трассе вычислений, чтобы убедиться в её корректности.
Доказатель отвечает на запросы, отправляя следующие данные:
1. Значения полинома $P(x)$ в запрошенных точках.
2. Аутентификационные пути (Merkle-доказательства) для подтверждения значений.
3. Сведения о соседних слоях FRI для проверки ограничения степени.

**Важно**:
- Верификатор проверяет только малую часть данных, но этого достаточно, чтобы убедиться с высокой вероятностью, что доказательство корректно.
- Безопасность основана на криптографически сильной случайности в запросах.


Эти этапы обеспечивают выполнение главных свойств ZK-STARKs: прозрачность, квантовую стойкость и масштабируемость. Следующие главы будут детализировать каждый из этих процессов с примерами и схемами.

---

## Заключение

ZK-STARKs представляют собой мощный инструмент для создания доказательств с нулевым знанием, который совмещает в себе прозрачность, масштабируемость и безопасность. В отличие от других подходов, таких как SNARKs, они полностью избавлены от необходимости в доверенной настройке, что делает их особенно подходящими для приложений, требующих высокого уровня доверия и защиты данных.

Основные преимущества ZK-STARKs:
1. **Прозрачность**: Доказательства не зависят от секретных параметров, а базируются на публично доступных криптографических примитивах (например, хеш-функциях).
2. **Квантовая устойчивость**: В отличие от алгоритмов, основанных на эллиптических кривых, ZK-STARKs устойчивы к атакам, использующим квантовые вычисления.
3. **Масштабируемость**: Протокол эффективно обрабатывает даже огромные вычислительные задачи, позволяя генерацию доказательств за приемлемое время, а проверка доказательств остаётся быстрой.
4. **Универсальность**: Подход применим к широкому спектру задач, от блокчейнов до приватных вычислений в распределённых системах.

Системы доказательств с нулевым знанием, такие как ZK-STARKs, находят применение в следующих областях:
- **Блокчейн**: Проверка транзакций без раскрытия входных данных, конфиденциальные смарт-контракты.
- **Обеспечение приватности**: Подтверждение подлинности данных без раскрытия их содержания.
- **Верификация вычислений**: Убедительное доказательство корректности вычислений для облачных сервисов или распределённых систем.

В заключение, ZK-STARKs делают ещё один шаг к созданию полностью децентрализованных, безопасных и конфиденциальных систем. С их помощью можно создать цифровой мир, где частная информация надёжно защищена, а доверие обеспечивается математической строгостью, а не зависимостью от третьих сторон.

Дальнейшие исследования и улучшения протоколов ZK-STARKs откроют ещё больше возможностей для интеграции нулевого знания в повседневные технологии, ускоряя переход к новому уровню приватности и безопасности.

