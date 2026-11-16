"""Build a small, varied set of English test passages (short / medium / long).

The grading passages are unseen, so we test on several styles. Run once:
    python texts/make_texts.py
"""
import os, re

HERE = os.path.dirname(os.path.abspath(__file__))

NARRATIVE = """The lighthouse keeper had not spoken to another person in eleven days when the boat appeared on the horizon. At first she took it for a trick of the light, one of the grey shapes the fog liked to make in the early morning, but it held its course and grew larger, and by nine o'clock she could make out a figure standing in the bow with one hand raised against the glare. She put the kettle on, then took it off again, then put it back on. It had been so long since she had needed to decide anything that the small choice felt enormous.

The figure turned out to be a boy of perhaps fifteen, thin and sunburnt, wearing a coat several sizes too large. He tied the boat to the iron ring at the foot of the steps with a knot that would not have held a rowboat in a bathtub, and she went down to retie it before she said a word to him. "You'll lose her like that," she said. "The tide comes in fast here, and it comes in angry."

"I didn't know where else to go," he said. He did not explain, and she did not ask. She had learned, in her years on the rock, that most questions answered themselves if you were patient enough to wait for them. She gave him tea and bread and the last of the marmalade, and he ate as though he had forgotten how, slowly at first and then all at once.

That night the wind rose, as it always did in the third week of October, and the two of them sat in the lamp room while the great lens turned above them, throwing its beam out across the black water. The boy watched it go round and round. "Does it ever stop?" he asked. "Not while I'm here," she said, and something in the way she said it made him smile for the first time.

In the morning the sea was flat and silver, and the boy was asleep in the chair with his coat pulled over his head. She let him sleep. There was work to do, there was always work to do, the brass to polish and the log to write and the wicks to trim, but for once she did it quietly, moving around him the way you might move around a sleeping animal you did not want to startle away."""

TECHNICAL = """A hash table stores key-value pairs in an array of buckets, using a hash function to compute the index at which each key should live. When two keys hash to the same bucket, a collision occurs, and the table must resolve it, either by chaining (keeping a small list in each bucket) or by open addressing (probing for the next free slot). The expected cost of a lookup depends on the load factor, the ratio of stored items to buckets. With a good hash function and a load factor below about 0.7, lookups, insertions and deletions all take constant time on average.

Resizing is the price of this speed. When the load factor crosses a threshold, the table allocates a larger array, typically twice the size, and reinserts every item. A single resize costs time proportional to the number of items, but because the size doubles each time, the total cost of all resizes over n insertions is still O(n), so the amortized cost per insertion remains constant. This argument, sometimes called the doubling trick, appears throughout the design of dynamic data structures.

Choosing a hash function is subtle. A function that is fast but poorly mixed will send structured keys, such as consecutive integers or strings sharing a long prefix, into a handful of buckets, and performance will collapse toward that of a linked list. Modern implementations often add a random seed to the hash at start-up, which defends against attackers who would otherwise craft inputs that all collide. Python, for example, randomizes string hashing by default for exactly this reason.

There are many refinements. Robin Hood hashing reduces the variance of probe lengths by letting an item that has travelled far from its home bucket displace one that has not. Cuckoo hashing gives each key two possible locations and guarantees worst-case constant-time lookups. Each trades some complexity in insertion for predictability in search, which matters most in systems where a single slow operation can stall an entire pipeline."""

DIALOGUE = """"Are you sure this is the right street?" asked Maya, peering at the map on her phone.
"Of course I'm sure," said Theo. "I've been here before. Twice."
"When?"
"Last summer. And maybe the summer before that. It looked different then."
"Different how?"
"There was a bakery on the corner. And fewer pigeons."
Maya looked up. There was no bakery on the corner. There was a laundromat, a closed bank, and a very large number of pigeons. "Theo, the restaurant closes in twenty minutes."
"Then we had better walk faster, hadn't we?"
They walked faster. They turned left, then right, then left again, and found themselves in a small square with a fountain that did not work and a statue of a man nobody had heard of. Theo stopped and frowned at the statue as though it had personally betrayed him.
"I don't remember him," he said.
"I'm calling a taxi."
"Wait, wait. Give me one more minute. I think it's just past the church."
"What church?"
Theo pointed. Behind the statue, half hidden by scaffolding, stood a small stone church, and beside it, glowing warmly in the dusk, a restaurant with a hand-painted sign. Maya stared at it, then at him.
"Twice," Theo said, grinning. "I told you."""

NEWS = """City officials announced on Tuesday that the long-delayed renovation of the downtown library will begin next spring, ending more than six years of debate over the building's future. The project, expected to cost $48 million, will add a children's wing, expand public computer access, and restore the original 1912 reading room, whose ceiling has been closed to visitors since a leak was discovered in 2019.

"This building belongs to everyone in the city," said Mayor Ellen Park at a press conference on the library steps. "For too long we treated it as a problem to be managed. Now we are treating it as an investment." The council approved the funding in a 7-2 vote last month after a series of public meetings that drew hundreds of residents.

Critics argued that the money would be better spent on neighborhood branches, several of which have reduced their hours because of staffing shortages. Councilmember David Ortiz, who voted against the plan, said he supported the renovation in principle but worried about the timeline. "We have heard these promises before," he said. "I want to see a shovel in the ground before I celebrate."

During construction, the main collection will move to a temporary location on Fifth Avenue, and the library expects to keep most services running. Officials said the work should be finished by the fall of 2029."""

POETRY = """The tide goes out, the tide comes in,
and leaves the sand as smooth as skin;
the gulls go up, the gulls come down,
and circle slow above the town.

I walked the shore when I was young
with salt and sunlight on my tongue,
and every shell I carried home
I thought was mine, and mine alone.

Now I am older, and I see
the sea was never kept by me;
it lent me what it chose to lend,
and took it back again, my friend."""


def license_excerpt(n):
    for p in ("/opt/anaconda3/lib/python3.12/LICENSE.txt",):
        if os.path.exists(p):
            t = open(p, encoding="utf-8").read()
            t = re.sub(r"[ \t]+", " ", t)
            return t[2000:2000 + n]
    return ""


def main():
    texts = {
        "short_poetry": POETRY[:300],
        "short_dialogue": DIALOGUE[:250],
        "short_news": NEWS[:400],
        "medium_narrative": NARRATIVE[:1200],
        "medium_technical": TECHNICAL[:1000],
        "medium_dialogue": DIALOGUE,
        "long_narrative": NARRATIVE,
        "long_technical": TECHNICAL,
        "long_news_mix": NEWS + "\n\n" + POETRY + "\n\n" + TECHNICAL[:1500],
        "long_legal": license_excerpt(4000),
    }
    for k, v in texts.items():
        if v:
            open(os.path.join(HERE, k + ".txt"), "w", encoding="utf-8").write(v)
            print(f"{k:20s} {len(v):6d} chars")


if __name__ == "__main__":
    main()
