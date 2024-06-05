import requests

def extract_barcode_positions_and_length(url):
    barcode_positions = []
    barcode_length = 0

    # Send a GET request to the website
    response = requests.get(url)

    # Split the response content into lines
    lines = response.text.split('\n')

    for i, line in enumerate(lines, start=1):
        if i == 5:
            # Extract the barcode length from line 5
            match = re.search(r'>C\((\d+),', line)
            if match:
                barcode_length = int(match.group(1))
        elif i >= 6 and not line.startswith('</pre>'):
            barcode_positions.append(line.strip())
        elif line.startswith('</pre>'):
            break
        print (barcode_positions, barcode_length)

    return barcode_positions, barcode_length