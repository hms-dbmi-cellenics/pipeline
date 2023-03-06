import requests

# Collect instance data from within instance
# https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/instancedata-data-retrieval.html
# List of available endpoints: https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/instancedata-data-categories.html

def get_instance_metadata(endpoints):

    print("Getting instance metadata...")

    result = {}

    INSTANCE_METADATA_ENDPOINT = "http://169.254.169.254/latest/meta-data/"

    def _query(metadata):
        url = f"{INSTANCE_METADATA_ENDPOINT}/{metadata}"
        response = requests.get(url)
        return response.text

    for endpoint in endpoints:
        result[endpoint] = _query(endpoint)

    return result