// Copyright 2016 Khaled Y. M. (https://github.com/kmehrunes)
/*
This file is part of GLex Genome Compression.

GLex is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License version 3
as published by the Free Software Foundation.

GLex is distributed in the hope that it will be useful, but 
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with GLex.  If not, see https://www.gnu.org/licenses/gpl.txt.
*/

using System;
using System.Collections.Generic;

namespace GLex
{
	// the possible encoding modes supported by GLex
	public enum EncodingMode: int {
		UNSIGNED_16_BITS = 0,
		UNSIGNED_32_BITS,
		UNSIGNED_64_BITS,
	}

	public enum StrandType {
		DNA,
		RNA
	}

	/// <summary>
	/// GLex encoding utility class, currently only supports DNA with raw A,C,G,T strings.
	/// </summary>
	public class GLexEncoding
	{
		static readonly Dictionary<char, byte> DNAEncodings = new Dictionary<char, byte> () {
			{ 'A', 0 },
			{ 'C', 1 },
			{ 'G', 2 },
			{ 'T', 3 }
		};

		static readonly Dictionary<byte, char> DNADecodings = new Dictionary<byte, char> () {
			{ 0, 'A' },
			{ 1, 'C' },
			{ 2, 'G' },
			{ 3, 'T' }
		};

		static readonly Dictionary<char, byte> RNAEncodings = new Dictionary<char, byte> () {
			{ 'A', 0 },
			{ 'C', 1 },
			{ 'G', 2 },
			{ 'U', 3 }
		};

		static readonly Dictionary<byte, char> RNADecodings = new Dictionary<byte, char> () {
			{ 0, 'A' },
			{ 1, 'C' },
			{ 2, 'G' },
			{ 3, 'U' }
		};

		public static byte EncodeSymbol (char nucleotide, StrandType strand = StrandType.DNA)
		{
			try {
				return strand == StrandType.DNA? DNAEncodings [nucleotide] : RNAEncodings [nucleotide];
			} catch {
				return 5; // used to mark failure
			}
		}

		public static char DecodeSymbol (byte encoding, StrandType strand = StrandType.DNA)
		{
			try {
				return strand == StrandType.DNA? DNADecodings [encoding] : RNADecodings [encoding];
			} catch {
				return (char)0; // used to mark failure
			}
		}

		/// <summary>
		/// Gets the maximum sequence length a certain encoding mode can represent.
		/// </summary>
		/// <returns>The sequence length.</returns>
		/// <param name="mode">The encoding mode</param>
		public static int MaxSequenceLength (EncodingMode mode)
		{
			ulong maxEncoding = 0;

			switch (mode) {
			case EncodingMode.UNSIGNED_16_BITS:
				maxEncoding = (ulong)ushort.MaxValue;
				break;

			case EncodingMode.UNSIGNED_32_BITS:
				maxEncoding = (ulong)uint.MaxValue;
				break;

			case EncodingMode.UNSIGNED_64_BITS:
				maxEncoding = ulong.MaxValue;
				break;
			}

			return LexicographicalDecoding (maxEncoding).Length;
		}

		/// <summary>
		/// Encodes a pattern lexicographically
		/// </summary>
		/// <returns>The lexicographical index of the pattern in a list of patterns of the same length</returns>
		/// <param name="pattern">Pattern.</param>
		public ulong LexicographicalEncoding (string pattern)
		{
			if (pattern == null) {
				throw new ArgumentException ("Pattern cannot be null");
			}

			if (pattern.Length == 0) {
				throw new ArgumentException ("Pattern cannot be empty");
			}

			return LexicographicalEncoding (pattern, pattern.Length);
		}

		/// <summary>
		/// Recursively find a lexicographical order of a pattern, which represents the index
		/// of that pattern assuming that all previous patterns are generated.
		/// </summary>
		/// <returns>The lexicographical index of the pattern in a list of patterns of length 'len'.</returns>
		/// <param name="pattern">Pattern.</param>
		/// <param name="len">The length of the pattern (used for recursion).</param>
		public static ulong LexicographicalEncoding (string pattern, int len, StrandType strand = StrandType.DNA)
		{
			byte encoding = EncodeSymbol (pattern [len - 1], strand);
			if (encoding > 4) {
				throw new ArgumentException ("The pattern is not a valid genome sequence");
			}

			if (len == 1) {
				return (ulong)encoding;
			}

			return 4 * LexicographicalEncoding (pattern, len - 1, strand) + encoding;
		}

		/// <summary>
		/// Decodes a lexicographical encoding.
		/// </summary>
		/// <returns>A pattern of length 'len' which has the index 'encoding' in a list of all pattern of length 'len'.</returns>
		/// <param name="encoding">The lexicographical encoding of the pattern (the index of the pattern).</param>
		/// <param name="len">The length of the pattern.</param>
		public static string LexicographicalDecoding (ulong encoding, int len, StrandType strand = StrandType.DNA)
		{
			if (len < 1 || encoding < 0) {
				throw new ArgumentException ("The encoding cannot be less than 0 and the length cannot be less than 1");
			}

			char[] decoded = new char[len];
			int index = len;
			ulong quotient = encoding, remainder;

			while (quotient > 0 || index > 0) {
				remainder = quotient % 4;
				quotient /= 4;
				decoded [--index] = DecodeSymbol ((byte)remainder, strand); // no need to check if the decoding failed, the remainder of dividing by 4 is always valid
			}

			return new string (decoded);
		}

		/// <summary>
		/// Decodes a lexicographical encoding.
		/// </summary>
		/// <returns>The shortest pattern which has the index 'encoding'. Note: All leading A's will be discarded</returns>
		/// <param name="encoding">The lexicographical encoding of the pattern (the index of the pattern).</param>
		public static string LexicographicalDecoding (ulong encoding, StrandType strand = StrandType.DNA)
		{
			if (encoding < 0) {
				throw new ArgumentException ("The encoding cannot be less than 0");
			}

			List<char> decoded = new List<char> ();
			ulong quotient = encoding, remainder;

			while (quotient > 0) {
				remainder = quotient % 4;
				quotient /= 4;
				decoded.Add (DecodeSymbol ((byte)remainder, strand)); // no need to check if the decoding failed, the remainder of dividing by 4 is always valid
			}

			return new string (decoded.ToArray());
		}
	}
}

